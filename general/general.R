# options(warn=0)

kDataDir <- "data/processed/"
kHomeDir   <- "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/"

library(data.table)
library(beeswarm)
source("general/Lists.R")


## TODO: Get ride of kDataDir and kHomeDir:
# Defines essential lists used throughout scripts:
# mirnas.all <- c("miR-1", "let-7a", "miR-155", "miR-124", "lsy-6")
# centered.sites <- c("11mer-m3.13", "12mer-m3.14", "11mer-m4.14", "12mer-m4.15")
# canonical.sites <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8")
kKmerSiteLists <- c("8mers", "9mers", "10mers", "11mers", "12mers")



kAgoStock <- data.frame(fread("SolveForKds/k_c_stockago.txt", sep="\t"), row.names=1)

kAgoRad <- data.frame(fread("SolveForKds/equilibrium_radioactivity_quants.txt", sep="\t"), row.names=1)
colnames(kAgoRad) <- gsub("^(.*)\\.(.*)$", colnames(kAgoRad), replace="\\1-\\2",
                          perl=TRUE)
colnames(kAgoRad)[ncol(kAgoRad)] <- "miR-7-23nt"

print(kAgoStock)
print(kAgoRad)


kAgoStockFunc <- function(mirna, logtrans=FALSE) {
  # This function fits the data to a table of radioactivity data, that gives the 
  # fraction bound of a 100 nM library during the RBNS reaction. This can be
  # used as an alternative way to assess the concentration of Ago, which is
  # provided in the original table from following the radioactivity from the
  # initial purification with radioactive guide strand.
  x_vals <- as.numeric(rownames(kAgoRad))
  y_vals <- kAgoRad[, mirna]*100 # This is done because the library is 100 nM
  if (logtrans) logstr <- "xy"
  else          logstr <- ""
  FitAgoData <- function(pars, logtrans=FALSE) {
    pars <- exp(pars)
    A <- pars[1]
    b <- pars[2]
    y_model <- A*x_vals + b
    if (logtrans) {
      y_vals <- log(y_vals)
      y_model <- log(y_model)
    }
    sum((y_vals - y_model)^2)
  }
  out <- exp(optim(c(1, 1), FitAgoData, logtrans=logtrans)$par)
}

kAgoStockNew <- sapply(kMirnasEquil, kAgoStockFunc)
kAgoStock[kMirnasEquil,"equilibrium"] <- kAgoStockNew[1, ]


## FUNCTIONS ###################################################################


## This is used almost everywhere throughout the R scripts as a convenience
## function to not have to re-assign the arguments when calling a function
## within another function, that share many of the same arguments.
SubfunctionCall <- function(f, ...) {
  # Check for the supplied arguments in "...":
  f.args.supplied <- list(...)
  # Remove these arguments from the list of all the arguments of function "f":
  f.arg.names <- setdiff(names(formals(f)), names(f.args.supplied))
  # Check for these arguments in the enviromnent:
  f.args <- lapply(f.arg.names, function(arg) {
    arg.temp <- try(dynGet(arg), silent=TRUE)
    if (class(arg.temp) != "try-error" & class(arg.temp) != "error") arg.temp
    else "remove"
  })
  names(f.args) <- f.arg.names
  # Combine the found arguments with the supplied and perform the function:
  f.args.list <- as.list(c(f.args[f.args != "remove"], f.args.supplied))
  do.call(f, f.args.list)
}


## Functions that load various types of data ###################################
EnsureDirectory <- function(dir_path) {
  if (!file.exists(dir_path)) {
    dir.create(dir_path)
  }
}


GetAnalysisPath <- function(mirna, experiment, condition, analysis_type, ext="",
                            suffix="txt") {
  dir_path <- sprintf("%s%s/%s/%s", kDataDir, mirna, experiment, analysis_type)
  print(dir_path)
  # dir_path <- sprintf("%s%s", kDataDir, sub_dir)
  EnsureDirectory(dir_path)
  sprintf("%s/%s%s.%s", dir_path, condition, ext, suffix)
}


GetMirnaCountData <- function(mirna, condition, unique=FALSE, no_adapter=FALSE,
                              no_marker=FALSE) {
  ext <- ""
  if (unique)  ext <- sprintf("%s_unique", ext)
  if (no_adapter) ext <- sprintf("%s_noadapter", ext)
  if (no_marker) ext <- sprintf("%s_nomarker", ext)
  path <- GetAnalysisPath(mirna, "AGO_purity", condition,
                           "AGO_pur_counts", ext=ext)
  file_counts <- read.table(path, sep="\t", header=TRUE, row.names=1)
  file_counts
}



SitesXCounts <- function(mirna, experiment="equilibrium", n_constant=5,
                         sitelist="resubmissionfinal", mirna.start=NULL,
                         multisite=FALSE, split16=NULL, uniq=FALSE,
                         start_mm=FALSE, stop_mm=FALSE, buffer=FALSE,
                         comp=FALSE, sc_k=6, m_k=6, new=FALSE, new2=FALSE,
                         include_zeros=FALSE, zeros_override=FALSE) {
  prog_sitelists <- c("programmed", "programmed_collapsed",
                      "programmed_suppcomp")
  ext <- sprintf("_%s", n_constant)
  print(sitelist)
  # if (sitelist %in% prog_sitelists & experiment %in% kExpsThreeP) {
  #   analysis_type <- "programmed_site_count_tables"
  #   if (class(start_mm) == "numeric" & class(stop_mm) == "numeric") {
  #     ext <- sprintf("%s_m%s.%smmsd", ext, start_mm, stop_mm)            
  #   }
  #   if (sitelist == "programmed_collapsed") {
  #     ext <- sprintf("%s_collapsed", ext)
  #   } else if (sitelist == "programmed_suppcomp") {
  #     ext <- sprintf("%s_suppcomp", ext)
  #   }
  if (multisite) {
    analysis_type <- "multisite_count_tables"
    ext <- sprintf("%s_%s", ext, sitelist)
  } else {
    analysis_type <- "site_count_tables"
    ext <- sprintf("%s_%s", ext, sitelist)
  }
  if (class(start_mm) == "numeric" & class(stop_mm) == "numeric") {
    ext <- sprintf("%s_m%s.%smmsd", ext, start_mm, stop_mm)            
  }
  condition <- "all_sites"
  if (sitelist %in% kKmerSiteLists) {
    ext <- sprintf("%s_%s-%s", ext, mirna.start, as.integer(mirna.start) + 3)
  }
  if (sitelist == "16mers") {
    ext <- sprintf("%s_%s", ext, split16)
  }
  if (uniq) {
    ext <- sprintf("%s_uniq", ext)
  }
  if (comp) {
    m_k <- min(m_k, c_k)
    ext <- sprintf("%s_comp_c%sm%s", ext, c_k, m_k)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  if (new) {
    ext <- sprintf("%s_new", ext)
  }
  if (new2) {
   ext <- sprintf("%s_new2", ext)
  }
  if (include_zeros) {
   ext <- sprintf("%s_zeros", ext)
  }
  file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  print(file)
  sXc <- data.frame(fread(file, sep="\t", header=TRUE, stringsAsFactors=FALSE),
                    row.names=1)
  # Remove "Seq" column, if it exists:
  sXc <- sXc[, grep("Seq", colnames(sXc), inv=TRUE)]
  # Remove the Leading `A` character from the AGO concentration samples.
  colnames(sXc) <- gsub("^A", colnames(sXc), replace="", perl=TRUE)
  # Replaces ".1" with ",1" at the end (this is really bad), because replicates
  # are named with ",1", ",2", ",3", etc.
  colnames(sXc) <- gsub("\\.([1234])$", colnames(sXc), replace=",\\1",
                        perl=TRUE)
  # This is a bandaid for the prior regex conversion because 0.4 is a standard
  # dilution used in the RBNS experiments.
  colnames(sXc)[which(colnames(sXc) == "0,4")] <- "0.4"
  # Check the two conditions, that the site exists in the input library:
  if (!(sitelist %in% c("12mers", "16mers")) & (!include_zeros) & (!zeros_override)) {
    sXc <- sXc[rowSums(sXc[, 3:ncol(sXc)], na.rm=TRUE) > 0,]
    if (!(multisite)) sXc <- sXc[sXc[, 2] > 0,]
  }
  # Conditional that orders the papercutoff3 sitelists based on the papercutoff
  # sitelist.
  if (sitelist == "papercutoff3") {
    sXc_alt <- SubfunctionCall(SitesXCounts, sitelist="papercutoff")
    sites <- rownames(sXc)
    sites_alt <- rownames(sXc_alt)
    sites <- sites_alt[which(sites_alt %in% sites)]
    sXc <- sXc[sites, ]
  }
  if (sitelist == "randthrp") {
    ind_none <- which(rownames(sXc) == "None")
    sXc <- rbind(sXc[-ind_none, ], sXc[ind_none, , drop=FALSE])
  }
  return(sXc)
}

EquilPars <- function(
  mirna, experiment="equilibrium", n_constant=5, sitelist="resubmissionfinal",
  uniq=FALSE, buffer=FALSE, mirna_start=FALSE, combined=TRUE, singleonly=FALSE,
  Xval=FALSE, L=FALSE, compcorrect=FALSE, wobble=FALSE, tpomit=FALSE,
  tpomit2=FALSE, tp2rep=FALSE, minkd=FALSE, AGOfixedbypass=FALSE, nbomitc=FALSE,
  sorted=FALSE, collapsemm=FALSE, sumseed=FALSE, start_mm=FALSE, stop_mm=FALSE,
  new=FALSE, new2=FALSE
) {
  analysis_type <- "kds_PAPER"
  condition <- ""
  if ((class(start_mm) == "numeric" | class(start_mm) == "integer") &
      (class(stop_mm) == "numeric" | class(stop_mm) == "integer")) {
      str_mm <- sprintf("m%s.%smmsd_", start_mm, stop_mm)
  } else {
    str_mm <- ""
  }
  ext <- sprintf("%s_%s%s", n_constant, str_mm, sitelist)
  if (sitelist %in% kKmerSiteLists) {
    ext <- sprintf("%s_%s-%s", ext, mirna_start, as.integer(mirna_start) + 3)
  }
  if (uniq) {
    ext <- sprintf("%s_uniq", ext)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  if (!combined) {
    ext <- sprintf("%s_nocombInput", ext)
  }
  if (singleonly) {
    ext <- sprintf("%s_singleonly", ext)
  }
  if (Xval) {
    ext <- sprintf("%s_Xval", ext)
  }
  if (L) {
    ext <- sprintf("%s_L%s", ext, L)
  }
  if (compcorrect) {
    ext <- sprintf("%s_compcorrect", ext)
    if (wobble) {
      ext <- sprintf("%s_wobble", ext)
    }
  }
  if (tpomit) {
    ext <- sprintf("%s_tpomit", ext)
  }
  if (tpomit2) {
    ext <- sprintf("%s_tpomit2", ext)
  }
  if (tp2rep) {
    ext <- sprintf("%s_rep%s", ext, tp2rep)
  }
  if (minkd) {
    minkd <- gsub("^(.*)e\\-0?(.*)$", sprintf("%.e", minkd),
                  replace="\\1e\\-\\2")
    ext <- sprintf("%s_minkd%s", ext, minkd)
  }
  if (AGOfixedbypass) {
    ext <- sprintf("%s_AGOfixedbypass", ext)
  }
  if (nbomitc) {
    ext <- sprintf("%s_nbomitc", ext)
  }
  if (collapsemm & sitelist %in% c("programmed_suppcomp", "bipartite_random_suppcomp")) {
    ext <- sprintf("%s_collapsemm", ext)
  }
  if (sumseed & sitelist %in% c("progthrp_suppcomp", "randthrp_suppcomp", "randthrp_comp")) {
    ext <- sprintf("%s_sumseed", ext)
  }
  if (new) {
    ext <- sprintf("%s_new", ext)
  }
  if (new2) {
    ext <- sprintf("%s_new2", ext)
  }
  

  # Add 'PAPER' to every file
  # TODO remove this extension from everything.
  ext <- sprintf("%s_PAPER", ext)
  # This is a suffix extension for the programmed libraries, looking at the
  # actual range of values fitted during the optimization to be able to compare
  # to the confidence intervals that I am calculating.
  if (sorted) {
    ext <- sprintf("%s_sorted", ext)
  }
  # Get path to file.
  params.file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  # system(sprintf("ls -l %s", params.file))
  # Get the parameter file.
  params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
                       sep= "\t", stringsAsFactors=FALSE))
  return(params)
}

## Mathematical functions ######################################################
Norm <- function(vector) {
  return(vector/sum(vector))
}


## String functions ############################################################
StrRev <- function(x) {
  # For a single string or a vector of strings, returns each element of the
  # string in reverse order.
  return(sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))
}

SeqComplement <- function(x, RNA=FALSE) {
  # For a single string or a vector of strings, returns the complementary RNA
  # sequence (not reversed!).
  complement.base <- c(U="A", T="A", G="C", C="G", N="N", V="B", B="V", D="H", H="D")
  if (RNA) complement <- c(complement.base, A="U")
  else complement <- c(complement.base, A="T")
  return(sapply(lapply(strsplit(x, NULL), function(y) complement[y]), paste,
                collapse=""))
}


RevComplement <- function(x, RNA=FALSE) {
  StrRev(SeqComplement(x, RNA=RNA))
}

ConvertTtoU <- function(site) {
  gsub("T", c(site), replace="U")
}


GetRevMirnaIndex <- function(mirna, ind) {
  return(nchar(kMirnaSeqs[mirna]) - ind + 1)
}


ConvertTtoUandMmtoX <- function(site) {
  gsub("mm", c(gsub("T", c(site), replace="U")), replace="x")
}


GetRemovedSites <- function(sXc) {
  # Identify unidentified sites ("CACACAC", "GCACTTTA", etc.)
  grep("^(?:A|C|G|T)+$", rownames(sXc), value=TRUE)
}


GetRepressionLinearModel <-  function(mirna, experiment="equilibrium",
                                      n_constant=5, sitelist="resubmissionfinal",
                                      buffer=FALSE, combined=TRUE, 
                                      bg_method=3, exrib=FALSE, new=FALSE, old=FALSE,
                                      intercept=FALSE, best_site=FALSE,
                                      single_site=FALSE, n_cutoff=-Inf,
                                      kd_cutoff=Inf) {
  repression_matrix <- SubfunctionCall(GetRepressionMatrix)
  # saved_names <- colnames(repression_matrix)
  # repression_matrix <- data.frame(repression_matrix, row.names=1)
  # colnames(repression_matrix) <- saved_names[-1]
  if (mirna == "miR-7-23nt") {
    mirna_trim <- "miR-7"
  } else {
    mirna_trim <- mirna
  }
  pars <- SubfunctionCall(EquilPars)
  if (single_site) {
    repression_matrix <- repression_matrix[which(rowSums(repression_matrix[,-ncol(repression_matrix)]) == 1),]
  }
  rownames(pars) <- gsub("_Kd", "", rownames(pars))
  sites <- grep(sprintf("(?:%s|None)", mirna_trim), rownames(pars), value=TRUE,
                invert=TRUE)
  sites <- colnames(repression_matrix[1:(ncol(repression_matrix) - 1)])
  site_occurence <- repression_matrix[,sites, drop=FALSE]

  n_sites <- colSums(site_occurence)

  kds <- pars[sites, ]$Mean

  sites <- sites[which((n_sites >= n_cutoff & kds <= kd_cutoff) | sites %in% kSeedSites)]
  n_sites <- n_sites[sites]

  if (length(sites) == 0) {
    return()
  }
  order_sites <- seq(length(sites))
  names(order_sites) <- sites
  site_occurence <- site_occurence[, sites]
  fc <- repression_matrix[["fc"]]/log(2)
  nosite_fc <- mean(fc[which(rowSums(site_occurence) == 0)])
  fc <- fc - nosite_fc
  if (intercept) {
    model <- lm(fc ~ ., data=site_occurence)    
  } else {
    model <- lm(fc ~ . - 1, data=site_occurence)

  }
  CIs <- confint(model)
  vals <-  model$fitted.values
  out <- summary(model)$coefficients[, 1:2, drop=FALSE]
  summary_global <<- summary(model)
  CIs <- CIs[which(n_sites > 0),]
  n_sites <- n_sites[which(n_sites > 0)]
  if (intercept) {
    out <- cbind(out[-1, , drop=FALSE], n_sites, CIs[-1, ])    
  } else {
    out <- cbind(out[, , drop=FALSE], n_sites, CIs)    
  }
  rownames(out) <- gsub("`", "", rownames(out))
  return(out)
}

GetRepressionMatrix <- function(mirna, sitelist, flanks=FALSE, bg_method=3, exrib=FALSE, new=FALSE, old=FALSE) {
  if (flanks) {
    flank_str = "_flank"
  } else {
    flank_str = ""
  }
  if (new) {
    new_str = "_new"
  } else if (old) {
    new_str = "_old"
  } else {
    new_str = ""
  }
  if (exrib) {
    exrib_str = "_exrib"
  } else {
    exrib_str = ""
  }
  file_path <- sprintf("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/repression_hela_cs/lin%s_model_df/%s_%s%s%s.txt",
                                     mirna, flank_str, sitelist, bg_method, exrib_str, new_str)
  print(file_path)
  repression_matrix <- fread(file_path, sep="\t")
  saved_names <- colnames(repression_matrix)
  repression_matrix <- data.frame(repression_matrix, row.names=1)
  colnames(repression_matrix) <- saved_names[-1]
  repression_matrix
}


################################################################################
## FUNCTIONS RELATING TO FLANKING DINULCEOTIDE DATA ############################
################################################################################

SiteFlanksXCounts <- function(mirna, site, experiment="equilibrium",
                              n_constant=5, sitelist="resubmissionfinal",
                              uniq=FALSE, buffer=FALSE) {
  analysis_type <- "site_count_tables"
  condition <- paste0(site, "_flanking")
  ext <- sprintf("_%s_%s", n_constant, sitelist)
  if (uniq) {
    ext <- spring("%s_uniq", ext)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  sitesXcounts <- read.table(file, stringsAsFactors=FALSE)
  colnames.temp <- colnames(sitesXcounts)[3 : ncol(sitesXcounts)]
  colnames.new <- sapply(colnames.temp, function(name) {
    return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
  })
  colnames(sitesXcounts)[3 : ncol(sitesXcounts)] <- colnames.new
  return(sitesXcounts)
}


SitesAndSingleFlanksXCounts <- function(mirna, site, experiment="equilibrium",
                                        n_constant=5,
                                        sitelist="resubmissionfinal",
                                        buffer=FALSE) {
  sXc <- SubfunctionCall(SitesXCounts)
  fXc <- SubfunctionCall(SiteFlanksXCounts)
  site.ind <- which(rownames(sXc) == site)
  left.flanks <- substr(rownames(fXc), 1, 2)
  right.flanks <- substr(rownames(fXc), 4, 5)
  rownames(fXc) <- paste0(site, "|", rownames(fXc))
  rbind(sXc[0:(site.ind - 1), ],
        fXc,
        sXc[(site.ind + 1): nrow(sXc), ])
}

GetFlankKds <- function(mirna, site, experiment="equilibrium", n_constant=5,
                        sitelist="resubmissionfinal", combined=TRUE,
                        buffer=FALSE) {

  analysis_type <- "kds_PAPER"
  condition <- ""
  ext <- sprintf("%s_%s_%s", n_constant, sitelist, site)
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  if (!combined) {
    ext <- sprintf("%s_nocombInput", ext)
  }
  ext <- sprintf("%s_PAPER", ext)
  params.file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  data <- fread(params.file, fill=TRUE,header=TRUE,
                       stringsAsFactors=FALSE, showProgress=FALSE)
  colnames(data) <- c("", colnames(data)[-ncol(data)])
  out <- data.frame(data, row.names=1)
  out
}


GetFullMirnaSiteFlankKds <- function(mirna, site, experiment="equilibrium",
                                     n_constant=5, sitelist="resubmissionfinal",
                                     buffer=FALSE, combined=TRUE) {
  tick <<- tick + 1
  out <- setNames(rep(NaN, length(kFlanks)), kFlanks)
  f_kds <- SubfunctionCall(GetFlankKds)
  site <- gsub("\\(", site, replace="\\\\(")
  site <- gsub("\\)", site, replace="\\\\)")
  f_kds <- f_kds[grep(sprintf("^%s\\|", site), rownames(f_kds), perl=TRUE), ]
  rownames(f_kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), replace="\\1",
                          rownames(f_kds), perl=TRUE)
  flanks <- rownames(f_kds)
  out[flanks] <- f_kds[flanks, ]$Mean
  return(out)
}

GetFullMirnaSiteFlankKds_CI <- function(mirna, site, experiment="equilibrium",
                                     n_constant=5, sitelist="resubmissionfinal",
                                     buffer=FALSE, combined=TRUE) {
  tick <<- tick + 1
  out <- setNames(rep(NaN, length(kFlanks)), kFlanks)
  out <- cbind(out, out)
  f_kds <- SubfunctionCall(GetFlankKds)
  site <- gsub("\\(", site, replace="\\\\(")
  site <- gsub("\\)", site, replace="\\\\)")
  f_kds <- f_kds[grep(sprintf("^%s\\|", site), rownames(f_kds), perl=TRUE), ]
  rownames(f_kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), replace="\\1",
                          rownames(f_kds), perl=TRUE)
  flanks <- rownames(f_kds)
  colnames(out) <- c("Upper_CI", "Lower_CI")
  out[flanks, 1] <- f_kds[flanks, ]$Upper_CI
  out[flanks, 2] <- f_kds[flanks, ]$Lower_CI 
  return(out)
}


GetFlankLinearModel <- function(experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                combined=TRUE, leaveout=FALSE) {
  
  if (leaveout != FALSE) {
    kMirnasEquil <- setdiff(kMirnasEquil, leaveout)
  }
  kd.data <- rbindlist(lapply(kMirnasEquil, function(mirna) {
    # if (mirna %in% c("miR-7-23nt")) {
    #   combined <- FALSE
    # }
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
      sitelist <- "resubmissionfinal6merm8"
    }
    if (mirna == "miR-1") {
      buffer <- TRUE
    }
    rbindlist(lapply(kSeedSites, function(site) {
      # Get each site and miRNA Kd table:
      kds <- SubfunctionCall(GetFlankKds)
      kds <- kds[grep(sprintf("^%s\\|", site), rownames(kds), perl=TRUE), ]
      rownames(kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), "\\1", rownames(kds), perl=TRUE)
      flank.5p <- ExtractSingleFlanks(rownames(kds))
      flank.3p <- ExtractSingleFlanks(rownames(kds), side="3p")
      data.frame(logkd=log(kds$Mean),
                 mirna=rep(mirna, nrow(kds)),
                 site=rep(site, nrow(kds)),
                 f5p1=substr(flank.5p, 1, 1),
                 f5p2=substr(flank.5p, 2, 2),
                 f3p1=substr(flank.3p, 1, 1),
                 f3p2=substr(flank.3p, 2, 2))
    }))
  }))
  kd_data_global <<- kd.data
  out <- lm(logkd ~ mirna*site + f5p1*f5p2 + f3p1*f3p2 - 1, data=kd.data)
  out
}

ExtractSingleFlanks <- function(flank.names, side="5p") {
  if (side == "5p") {
    replace.string <- "\\1"
  } else if (side == "3p") {
    replace.string <- "\\2"
  }
  return(gsub("^(.{2,2}).*(.{2,2})$", flank.names, replace=replace.string,
         perl=TRUE))
}

GetSingleMirnaModel <- function(mirna, experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                combined=TRUE) {
  lm(logkd ~ site - 1, data=SubfunctionCall(GetMirnaFlankDataFrame))
}


LeaveOneOutFlankModel <- function(mirna, experiment="equilibrium",
                                  n_constant=5, sitelist="resubmissionfinal",
                                  combined=TRUE, xpos=20, ypos=20) {
  # Get the flanking dinucleotide model trained on the other 5 miRNAs:
  flank.model <- SubfunctionCall(GetFlankLinearModel, leaveout=mirna)
  flank.model_temp <<- flank.model$model[1:256, ]
  flank.model <<- flank.model
  if (mirna == "miR-7-23nt") {
    experiment <- "equilibrium2_nb"
    sitelist <- "resubmissionfinal6merm8"
  }
  if (mirna == "miR-1") {
    buffer <- TRUE
  }
  # Get the site coefficients for this miRNA:
  site.mir.model <- SubfunctionCall(GetSingleMirnaModel)
  site.mir.model <<- site.mir.model
  # Get the site coefficients:
  site.coeffs <- site.mir.model$coefficients
  names(site.coeffs) <- gsub("^site(.*)$", names(site.coeffs),
                             replace="\\1", perl=TRUE)
  leave_out_coefficients <- c(site.mir.model$coefficients,
                                   flank.model$coefficients[grep("^f(?:5|3p)",
                                                                 names(flank.model$coefficients),
                                                                 perl=TRUE, value=TRUE)])


  switch_coeffs <- flank.model$coefficients[grep("^f(?:5|3p)",
                                                                 names(flank.model$coefficients),
                                                                 perl=TRUE, value=TRUE)]

  data.flanks <- SubfunctionCall(GetMirnaFlankDataFrame)
  lm_site <- lm(logkd ~ site - 1, data=data.flanks)
  lm_new <- lm(logkd ~ site + f5p1*f5p2 + f3p1*f3p2 - 1, data=data.flanks)
  lm_new$coefficients <- leave_out_coefficients
  lm_flanks_only_new <- lm(logkd ~ f5p1*f5p2 + f3p1*f3p2, data=data.flanks)
  lm_flanks_only_switch <- lm_flanks_only_new
  lm_flanks_only_switch$coefficients <- c(lm_flanks_only_new$coefficients[1], switch_coeffs)
  names(lm_flanks_only_switch$coefficients)[1] <- "(Intercept)"
  flanks_only_pred <<- predict(lm_flanks_only_new, data=data.flanks)
  flanks_only_switch_pred <<- predict(lm_flanks_only_switch, data=data.flanks)
  log_kds_new <- data.flanks$logkd - flanks_only_switch_pred
  log_site_new <- data.flanks$logkd - flanks_only_pred
  data.flanks_alt <- data.flanks
  data.flanks_alt$logkd <- log_kds_new
  lm_site_corrected <- lm(logkd ~ site - 1, data=data.flanks_alt)

  x_vars <- data.flanks$logkd

  y_vars <- predict(lm_new, data=data.flanks)


  lm_new$coefficients <- c(lm_site_corrected$coefficients,
                           lm_flanks_only_switch$coefficients[-1])

  y_vars_new <- predict(lm_new, data=data.flanks)

  # Get the flank coefficients:
  flank.coeffs <- SubfunctionCall(ParseFlankModelCoefs, model=flank.model)
  flank.lin <- flank.coeffs[[1]]
  flank.5int <- flank.coeffs[[2]]
  flank.3int <- flank.coeffs[[3]]
  data.flanks <- SubfunctionCall(GetMirnaFlankDataFrame)
  data.flanks <<- data.flanks
  test.model <- apply(data.flanks, 1, function(row) {
    site <- paste0("site", row["site"])
    f5p1 <- row["f5p1"]
    f5p2 <- row["f5p2"]
    f3p1 <- row["f3p1"]
    f3p2 <- row["f3p2"]
    site.val <- site.coeffs[row["site"]]
    f5p1.val <- flank.lin[f5p1, "f5p1"]
    f5p2.val <- flank.lin[f5p2, "f5p2"]
    f3p1.val <- flank.lin[f3p1, "f3p1"]
    f3p2.val <- flank.lin[f3p2, "f3p2"]
    if (f5p1 != "A" & f5p2 != "A") {
      int5.val <- flank.5int[paste0(2, f5p2), paste0(1, f5p1)]    
    } else {
      int5.val <- 0
    }
    if (f3p1 != "A" & f3p2 != "A") {
      int3.val <- flank.3int[paste0(2, f3p2), paste0(1, f3p1)]    
    } else {
      int3.val <- 0
    }
    sum(site.val, f5p1.val, f5p2.val, f3p1.val, f3p2.val, int5.val, int3.val)
  })
  alt_out <- test.model - mean(test.model, na.rm=TRUE) + mean(data.flanks$logkd, na.rm=TRUE)
  y_vars <- y_vars - mean(y_vars, na.rm=TRUE) + mean(data.flanks$logkd, na.rm=TRUE)
  colors <- kSiteColors[as.character(data.flanks$site)]

  y_vars_new
}


GetMirnaFlankDataFrame <- function(mirna, experiment="equilibrium",
                                   n_constant=5, sitelist="resubmissionfinal",
                                   buffer=FALSE, combined=TRUE) {
  rbindlist(lapply(kSeedSites, function(site) {
      # Get each site and miRNA Kd table:
      kds <- SubfunctionCall(GetFlankKds)
      kds <- kds[grep(sprintf("^%s\\|", site), rownames(kds), perl=TRUE), ]
      rownames(kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), "\\1", rownames(kds), perl=TRUE)
      flank.5p <- ExtractSingleFlanks(rownames(kds))
      flank.3p <- ExtractSingleFlanks(rownames(kds), side="3p")
      data.frame(logkd=log(kds$Mean),
                 mirna=rep(mirna, nrow(kds)),
                 site=rep(site, nrow(kds)),
                 f5p1=sapply(flank.5p, substr, start=1, stop=1),
                 f5p2=sapply(flank.5p, substr, start=2, stop=2),
                 f3p1=sapply(flank.3p, substr, start=1, stop=1),
                 f3p2=sapply(flank.3p, substr, start=2, stop=2))
  }))
}


ParseFlankModelCoefs <- function(model) {
  coefs <- c(f5p1A=0, f5p2A=0, f3p1A=0, f3p2A=0, model$coefficients)
  coefs_error <- c(f5p1A=0, f5p2A=0, f3p1A=0, f3p2A=0, summary(model)[[4]][,2])
  CI_lower <- c(f5p1A=0, f5p2A=0, f3p1A=0, f3p2A=0, confint(model)[, 1])
  CI_upper <- c(f5p1A=0, f5p2A=0, f3p1A=0, f3p2A=0, confint(model)[, 2])
  lin <- sapply(c("f5p1", "f5p2", "f3p1", "f3p2"), function(string) {
    grep.key <- paste0("^", string, ".{1}$")
    out <- coefs[grep(grep.key, names(coefs), perl=TRUE)]
    names(out) <- sapply(names(out), substr, start=5, stop=5)
    out
  })

  lin_error_upper <- sapply(c("f5p1", "f5p2", "f3p1", "f3p2"), function(string) {
    grep.key <- paste0("^", string, ".{1}$")
    out <- CI_upper[grep(grep.key, names(CI_upper), perl=TRUE)]
    names(out) <- sapply(names(out), substr, start=5, stop=5)
    out
  })


  lin_error_lower <- sapply(c("f5p1", "f5p2", "f3p1", "f3p2"), function(string) {
    grep.key <- paste0("^", string, ".{1}$")
    out <- CI_lower[grep(grep.key, names(CI_lower), perl=TRUE)]
    names(out) <- sapply(names(out), substr, start=5, stop=5)
    out
  })

  ints <- lapply(c(5, 3), function(flank){
    sapply(c("1C", "1G", "1T"), function(string) {
    grep.key <- paste0("^f", flank, "p", string, ":")
    out <- coefs[grep(grep.key, names(coefs), perl=TRUE)]
    len <- nchar(names(out))
    names(out) <- sapply(names(out), substr, start=len - 1, stop=len)
    out
    })  
  })
  list(lin=lin, int.5p=ints[[1]], int.3p=ints[[2]],
       lin_error_upper=lin_error_upper, lin_error_lower=lin_error_lower)
}


GetFlankLMCoefs <- function(experiment="equilibrium", n_constant=5,
                            sitelist="resubmissionfinal") {
  model <- SubfunctionCall(GetFlankLinearModel)
  SubfunctionCall(ParseFlankModelCoefs)
}


GetPairingFlankData <- function(mirna, site, condition,
                                experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                mir.start=1, mir.stop=14, noconstant=FALSE,
                                alt_mir_exp_cond=NULL,
                                old=FALSE) {
  ############# THIS FUNCTION TAKES DATA GENERATED BY ##########################
  # folder = "/structural_analysis_PAPER_realfinal/"
  # file = paste0(condition, "_", n_constant, "_", sitelist, "_", mir.start, "-", mir.stop,
  #               collapse="")
  # path <- paste0(kSolexaDir, mirna, "/", experiment, folder, site, "/", file,
  #                ".txt")
  # print(path)
  # if (old) {
  #   analysis_type <- "structural_analysis"  
  # } else {
  analysis_type <- "structural_analysis"  
  # }
  ext <- sprintf("_%s_%s_%s-%s", n_constant, sitelist, mir.start, mir.stop)
  if (!is.null(alt_mir_exp_cond)) {
    ext <- sprintf("%s_alt_%s", ext, alt_mir_exp_cond)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  if (noconstant == TRUE) {
    ext <- sprintf("%s_noconstant", ext)
  }
  path <- SubfunctionCall(GetAnalysisPath,
                              condition=sprintf("%s/%s", site, condition),
                              ext=ext)
  system(sprintf("ls -l %s", path))
  # data <- read.table(path, header=TRUE)
  data <- data.frame(read.table(path, header=TRUE))
  return(data)
}


# # Load relevant libraries:
# library(beeswarm)
# library(gplots)
# # library(wrswoR)
# library(deSolve)


# mirna_sequences <- c("UGGAAUGUAAAGAAGUAUGUAU",
#                      "UGAGGUAGUAGGUUGUAUAGUU",
#                      "UUAAUGCUAAUCGUGAUAGGGGU",
#                      "UAAGGCACGCGGUGAAUGCCAA",
#                      "UUUUGUAUGAGACGCAUUUCGA",
#                      "UGGAAGACUAGUGAUUUUGUUGUUU")
# names(mirna_sequences) <- c("miR-1",
#                             "let-7a",
#                             "miR-155",
#                             "miR-124",
#                             "lsy-6",
#                             "miR-7")
# RC_vector <- c("A", "C", "G", "T")
# names(RC_vector) <- c("U", "G", "C", "A")


# # Define colormaps for plotting:
# kEquilCondColors <- c("black", "black", "gray60", "red", "orange", "forestgreen",
#                       "blue", "violet")
# names(kEquilCondColors) <- c("I", "I_combined", "0", "0.4", "1.26", "4", "12.6",
#                              "40")

# kBaekColors <- c("purple1", "firebrick", "blue", "cyan", "purple3", "purple2",
#                  "lightblue", "darkslategrey", "darkslategray4",
#                  "darkslategray3", "darkslategray2", "black")
# names(kBaekColors) <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "7mer-m3.9",
#                         "6mer-m8", "6mer-A1", "CDNST 1", "CDNST 2", "CDNST 3",
#                         "CDNST 4", "None")

# kMirnaColors <- c("deepskyblue2", "black", "red", "forestgreen", "purple")
# names(kMirnaColors) <- mirnas.all

# kNucleotideColors <- c("blue", "green", "purple", "red")
# names(kNucleotideColors) <- c("A", "T", "C", "G")


# repression.df <- data.frame(read.table("RepressionData/best_flanking_kds_and_repression_3pseq.txt",header=TRUE,na.strings="",sep="\t", stringsAsFactors=FALSE))
# repression_old.df <- data.frame(read.table("RepressionData/all_flanking_kds_and_repression.txt",header=TRUE,na.strings="",sep="\t", stringsAsFactors=FALSE))
# repression_new <- repression.df
# repression2.df <- data.frame(read.table("RepressionData/all_flanking_kds_and_repression_3pseq.txt",header=TRUE,na.strings="",sep="\t", stringsAsFactors=FALSE))






# FilterSingleCanonicalSites <- function(repression_dataframe) {
#   apply(repression_dataframe, 1, function(row) {
#     print(row)
#     return(row)
#     })
# }
# tick <- 0


# GiveTargetSequence <- function(mirna, start, stop) {
#   sequence <- mirna_sequences[mirna]
#   sequence_window <- unlist(strsplit(sequence, split = ""))[start:stop]
#   sequence_converted <- RC_vector[rev(sequence_window)]
#   sequence_new <- paste0(sequence_converted, collapse = "")
#   return(sequence_new)
# }


# kSiteColors <- read.table(
#   "general/site_info.txt",
#   row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 1, drop = FALSE]
# names.temp <- rownames(kSiteColors)
# kSiteColors <- c(kSiteColors[,1])
# names(kSiteColors) <- names.temp

# kSiteCategoryColors <- read.table(
#   "general/site_info.txt",
#   row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)[, 6, drop = FALSE]
# names.temp <- rownames(kSiteCategoryColors)
# kSiteCategoryColors <- kSiteCategoryColors[,1]
# names(kSiteCategoryColors) <- names.temp


# stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")


# ConvertRColortoRGB <- function(color, alpha) {
#   rgb.vals <- c(col2rgb(color)) / 255
#   r <- rgb.vals[1]
#   g <- rgb.vals[2]
#   b <- rgb.vals[3]
#   return(rgb(r, g, b, alpha=alpha))
# }

# GetSingleFlankPosition <- function(flanks, position) {
#   if (position > 2) {
#     position <- position + 1
#   }
#   output <- sapply(flanks, function(flank_name) {
#     unlist(strsplit(flank_name, split=""))[position]
#     })
# return(output)

# }

# kPlotParameters <- list(
#   cex.main  = 1.5,
#   lwd       = 2,
#   pch       = 20,
#   cex.lab   = 1.5,
#   cex.axis  = 1.5,
#   ann       = FALSE,
#   font      = 1,
#   las       = 1,
#   mar       = c(5, 5, 4, 2) + 0.1,
#   font.main = 1,
#   bty       = "n",
#   mgp       = c(2.2, 1, 0))

# kPlotParameters <- list(
#   cex.main  = 1,
#   lwd       = 1.5,
#   lheight   = 1,
#   pch       = 19,
#   cex.lab   = 1,
#   lab       = c(5, 5, 2),
#   cex.axis  = 1,
#   ann       = FALSE,
#   font      = 1,
#   las       = 1,
#   mar       = c(3, 3, 2, 1),
#   font.main = 1,
#   tcl       = -0.2,
#   bty       = "n",
#   mgp       = c(0.5, 0.3, 0))




# Cumul <- function(vector) {
#   norm <- Norm(vector)
#   tot <- 0
#   out <- sapply(norm, function(x){
#     tot <<- tot + x
#     return(tot)

#     })
#   return(out)
# }

# Logistic <- function(vector, max) {
#   return(max/(1 + exp(-vector)))
# }

# Logit <- function(vector, max) {
#   return(-log(max/vector - 1))
# }


# ## FUNCTIONS USED IN EQUILIBRIUM MODELING (4 total)
# # 1.________
# GetOccupancy <- function(a, kds) {
#   # Generates occupancy matrix for the entire input matrix which has rows
#   # that are all possible flnaking nucleotide combinations, and columns that
#   # are the probabilities of being unpaired in the window.
#   #
#   # Args:
#   # a: The concentration of free AGO in the binding reaction
#   # kds: a list of all kds corresponding to the individual site types including
#   # Returns:
#   oc <- a/(a + kds)
#   return(oc)
# }

# GetOccupancyMulti <- function(a, kd.set) {
#   # Args:
#   # a: The concentration of free AGO in the binding reaction
#   # kds: A list of all
#   oc <- sapply(kd.set, function(kd) {
#     print(kd)
#     return(1 - 1/prod(1 + a/kd))
#   })
#   return(oc)
# }


# kd.set <- c(0.001, 0.002)

# data <- GetSitesXCounts

# ModelFunctionNew <- function(pars) {
#   # print(c.totals[1,1])
#   # print(c.I.tots[1])
#   # Split up the parameters into the kd and background parameters.
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


# GetAverageOfReplicates <- function(time, times, data) {

#   return(rowMeans(data[, which(times == time), drop=FALSE]))
# }



# GetContamResidual <- function(b, kd.b, l, A, a, B, exp=1){
#   # a: The concentration of free AGO in the binding reaction
#   # b: The concentration of free contaminant in the binding reaction
#   # A: The concentration of free AGO in the binding reaction
#   # B: The concentration of total contaminant in the binding reaction
#   # kd.b: The Kd value for the contaminant for all RNA.
#   res <- (b^2 + (sum(l) - A - B + a + kd.b)*b - B*kd.b)^exp
#   return(res)
# }

# GetFreeContam <- function(kds.a, kd.b, l, A, a, B) {
#   solution <- NaN
#   try(solution <- uniroot(GetContamResidual, c(0, B), kd.b=kd.b, l=l, A=A, B=B,
#                           a=a, tol=1e-6*.Machine$double.eps^0.25)$root,
#       silent=TRUE)
#   if (is.na(solution)) {
#     solution <- optimize(GetContamResidual, c(0, B), kd.b=kd.b, l=l, A=A, B=B,
#                          a=a, exp=2, tol=1e-6*.Machine$double.eps^0.25)$minimum
#   }
#   return(solution)
# }

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

# GetOccupanciesContaminant <- function(a, b, kds.a, kd.b) {
#   # Calculate common denomenator:
#   ocs.den <- kds.a*kd.b + kds.a*b + kd.b*a
#   return(list(a*kd.b/ocs.den, b*kds.a/ocs.den))
# }

# # 2._______
# GetFreeResidual <- function(a, kds, l, A, exp=1) {
# #
#   oc <- GetOccupancy(a, kds)
#   res <- (A - a - sum(oc*l))^exp
#   return(res)
# }

# GetFreeResidualMulti <- function(a, kd.set, l, A, exp=1) {
#   oc <- sapply(kd.set, function(kds){
#     return(sum(GetOccupancy(a, kds)))
#   })
#   res <- (A - a - sum(oc*l))^exp
#   return(res)
# }

# # 4._______



# GetFreeAgo <- function(kds, l, A, multi=FALSE) {
#   if (A > 0) {
#     if (multi == TRUE) {
#       ResidualFunction <- GetFreeResidualMulti
#     } else {
#       ResidualFunction <- GetFreeResidual
#     }
#     a <- NaN
#     try(a <- uniroot(ResidualFunction, c(0, A), kds=kds, l=l, A=A,
#                      tol=1e-4*.Machine$double.eps^0.25)$root)
#     if (is.na(a)) {
#       a <- optimize(ResidualFunction, c(0, A), kds=kds, l=l, A=A,
#                      tol=1e-5*.Machine$double.eps^0.25)$minimum
#     } 
#     return(a)
#   } else {
#     return(0)
#   }
# }

# GetBoundRNA <- function(kds, l, A, multi=FALSE) {
#   if (multi == TRUE) {
#     OccupancyFunction <- GetOccupancyMulti
#   } else {
#     OccupancyFunction <- GetOccupancy
#   }
#   a <- GetFreeAgo(kds, l, A, multi=multi)
#   x <- OccupancyFunction(a, kds)*l
#   return(x)
# }


# RemovePrefixFromString <- function(name, prefix) {
#   if (substr(name, 1, nchar(prefix)) == prefix) {
#     return(substr(name, nchar(prefix) + 1, nchar(name)))
#   } else {
#     return(name)
#   }
# }



# ## Functions that load data. (__ total) # 1.___________
# GetSitesXCounts <- function(mirna, exp, n_constant, sitelist, mirna.start=NULL,
#                             mirna.stop=NULL) {
#   file.prefix <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/", exp,
#                         "/site_count_tables/all_sites_", n_constant, "_",
#                         sitelist)
#   if (sitelist %in% kKmerSiteLists) {
#     file.suffix <- paste0("_", mirna.start, "-", mirna.stop, ".txt")
#     col.start <- 3
#   } else {
#     file.suffix <- ".txt"
#     col.start <- 4     
#   }
#   file <- paste0(file.prefix, file.suffix)
#   sXc <- read.table(file, stringsAsFactors=FALSE)
#   colnames(sXc) <- gsub("^A", colnames(sXc), replace = "")
#   return(sXc)
# }

# GetMultiSitesXCounts <- function(mirna, exp, n_constant, sitelist,
#                             mirna.start=NULL, mirna.stop=NULL) {
#   if (sitelist %in% kKmerSiteLists) {
#    sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                              mirna, "/", exp,
#                              "/multisite_count_tables/all_sites_",
#                              n_constant, "_", sitelist, "_", mirna.start, "-",
#                              mirna.stop, ".txt")
#     sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
#     colnames.temp <- colnames(sitesXcounts)[3 : ncol(sitesXcounts)]
#     colnames.new <- sapply(colnames.temp, function(name) {
#         return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
#       })
#     colnames(sitesXcounts)[3 : ncol(sitesXcounts)] <- colnames.new
#   } else {
#     sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                          mirna, "/",exp,"/multisite_count_tables/all_sites_",
#                          n_constant,"_", sitelist, ".txt")
#     sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
#     colnames.temp <- colnames(sitesXcounts)[4 : ncol(sitesXcounts)]
#     colnames.new <- sapply(colnames.temp, function(name) {
#         return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
#       })
#     colnames(sitesXcounts)[4 : ncol(sitesXcounts)] <- colnames.new
#   }
#   return(sitesXcounts)
# }


# GetSitesXCountsUnique <- function(mirna, experiment, start, stop, sitelist,
#                             mirna.start=NULL, mirna.stop=NULL) {
#   if (sitelist %in% kKmerSiteLists) {
#    sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                              mirna, "/", experiment,
#                              "/full_site_count_tables_unique/all_sites_",
#                              start, "-", stop, "_", sitelist, "_", mirna.start, "-",
#                              mirna.stop, ".txt")
#     sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
#     colnames.temp <- colnames(sitesXcounts)[2 : ncol(sitesXcounts)]
#     colnames.new <- sapply(colnames.temp, function(name) {
#         return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
#       })
#     colnames(sitesXcounts)[2 : ncol(sitesXcounts)] <- colnames.new
#   } else {
#     sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                          mirna, "/",experiment,"/full_site_count_tables_unique/all_sites_",
#                          start, "-", stop,"_", sitelist, ".txt")
#     sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
#     colnames.temp <- colnames(sitesXcounts)[3 : ncol(sitesXcounts)]
#     colnames.new <- sapply(colnames.temp, function(name) {
#         return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
#       })
#     colnames(sitesXcounts)[3 : ncol(sitesXcounts)] <- colnames.new
#   }
#   return(sitesXcounts)
# }



# GetSitesXCountsKinetics <- function(mirna, experiment, n_constant, sitelist) {
#   sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                        mirna, "/",experiment,"/site_count_tables/all_sites_",
#                        n_constant, "_", sitelist, "_pulse.txt")
#   sitesXcounts.pulse <- read.table(sites_file_name)
#   sitesXcounts.pulse <- sitesXcounts.pulse[,-4]
#   colnames(sitesXcounts.pulse)[3] <- "I_combined"
#   colnames(sitesXcounts.pulse) <- sapply(colnames(sitesXcounts.pulse),
#                                          RemovePrefixFromString, prefix="X")
#   colnames(sitesXcounts.pulse)[(
#   ncol(sitesXcounts.pulse)-2):ncol(sitesXcounts.pulse)] <- c(14400,"Equil-","0-")


#   sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                        mirna, "/",experiment,"/site_count_tables/all_sites_",
#                        n_constant, "_", sitelist, "_chase.txt")
#   sitesXcounts.chase <- read.table(sites_file_name)
#   sitesXcounts.chase <- sitesXcounts.chase[,-3]
#   colnames(sitesXcounts.chase)[3] <- "I_combined"
#   colnames(sitesXcounts.chase) <- sapply(colnames(sitesXcounts.chase),
#                                          RemovePrefixFromString, prefix="X")
#   colnames(sitesXcounts.chase)[(
#   ncol(sitesXcounts.chase)-2):ncol(sitesXcounts.chase)] <- c(14400,"Equil-","0-")

#   return(list(sitesXcounts.pulse,sitesXcounts.chase))
# }



# # 3.___________
# GetSiteFlanksXCounts <- function(mirna, experiment, n_constant,
#                                  sitelist, site) {
#   sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                        mirna, "/",experiment,"/site_count_tables/", site,
#                        "_flanking_", n_constant, "_", sitelist, ".txt")
#   print(sites_file_name)
#   sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
#   print(dim(sitesXcounts))
#   colnames.temp <- colnames(sitesXcounts)[3 : ncol(sitesXcounts)]
#   colnames.new <- sapply(colnames.temp, function(name) {
#     return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
#   })
#   colnames(sitesXcounts)[3 : ncol(sitesXcounts)] <- colnames.new
#   return(sitesXcounts)
# }

# GetSiteFlanksXCountsKinetics <- function(mirna, experiment, n_constant, sitelist, site) {
#   sites_file_name_pulse <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                        mirna, "/",experiment,"/site_count_tables/", site,
#                        "_flanking_pulse_", n_constant,"_", sitelist, ".txt")
#   sitesXcounts.pulse <- read.table(sites_file_name_pulse, stringsAsFactors=FALSE)
#   sitesXcounts.pulse <- sitesXcounts.pulse[,-3]

#   colnames(sitesXcounts.pulse)[2] <- "I_combined"
#   colnames(sitesXcounts.pulse) <- sapply(colnames(sitesXcounts.pulse),
#                                          RemovePrefixFromString, prefix="X")
#   colnames(sitesXcounts.pulse)[(
#   ncol(sitesXcounts.pulse)-2):ncol(sitesXcounts.pulse)] <- c(14400,"Equil-","0-")

#   sites_file_name_chase <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                        mirna, "/",experiment,"/site_count_tables/", site,
#                        "_flanking_chase_", n_constant, "_", sitelist, ".txt")
#   sitesXcounts.chase <- read.table(sites_file_name_chase, stringsAsFactors=FALSE)
#   sitesXcounts.chase <- sitesXcounts.chase[,-2]
#   colnames(sitesXcounts.chase)[2] <- "I_combined"
#   colnames(sitesXcounts.chase) <- sapply(colnames(sitesXcounts.chase),
#                                          RemovePrefixFromString, prefix="X")
#   colnames(sitesXcounts.chase)[(
#   ncol(sitesXcounts.chase)-2):ncol(sitesXcounts.chase)] <- c(14400,"Equil-","0-")

#   return(list(sitesXcounts.pulse, sitesXcounts.chase))
# }

# GetSiteKmersXCountsKinetics <- function(mirna, experiment, site, start, stop,
#                                  sitelist, k) {
#   sites_file_name_pulse <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                        mirna, "/",experiment,"/full_site_count_tables/", site,
#                        "_sitekmers_", start, "-", stop,"_", sitelist, "_k", k, "_pulse.txt")
#   sitesXcounts.pulse <- read.table(sites_file_name_pulse, stringsAsFactors=FALSE)
#   sites_file_name_chase <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#                        mirna, "/",experiment,"/full_site_count_tables/", site,
#                        "_sitekmers_", start, "-", stop,"_", sitelist, "_k", k, "_chase.txt")
#   sitesXcounts.chase <- read.table(sites_file_name_chase, stringsAsFactors=FALSE)

#   return(list(sitesXcounts.pulse, sitesXcounts.chase))
# }


# # 4.
# LinearModelInputFlanks <- function(flank_matrix,col) {
#   flanks <- matrix(sapply(c(1,2,4,5), function(index) {
#   sapply(rownames(flank_matrix), function(flank_name) {
#     unlist(strsplit(flank_name, split=""))[index]
#     })
#   }), nrow=nrow(flank_matrix), ncol=4, byrow=FALSE)
#   input_flanks <- data.frame(I = as.numeric(flank_matrix[,col]),
#                              f5p.o = flanks[, 1],
#                              f5p.i = flanks[, 2],
#                              f3p.i = flanks[, 3],
#                              f3p.o = flanks[, 4],
#                              stringsAsFactors = FALSE)
#   return(input_flanks)

# }



# ConvertTtoU <- function(site) {
#   sites.split <- strsplit(site, "")
#   out <- unlist(lapply(sites.split, function(site.split) {
#     ind.T <- grep("T",site.split)
#     site.split[ind.T] <- "U"
#     return(paste0(site.split, collapse = ""))
#     }))
#   return(out)
# }

# MakeIterationPlot <- function(out, type, extension="") {
#   setEPS()
#   postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                          "/figures/kds/", mirna, "/", type, "/iterations/",
#                          site, "_", k.c.stockago, extension, ".eps")
#             )
#   par(kPlotParameters)
#   x = seq(dim(out)[1])
#   probs = out[ ,"-logp"]
#   out.print <- out[ , seq(dim(out)[2] - 1)]
#   ys <- 10^c(floor(min(out.print)), ceiling(max(out.print)))

#   probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) +
#                   log10(ys[1])

#   out.print <- 10^out.print

#   plot(x , 10^probs.scaled, log='y', axes=FALSE, type="l", ylim=ys,
#        lwd=2, ann=FALSE,
#        col="black")
#   title(main=mirna, line=-1, adj=0.1)
#   title(main=site, col.main=kSiteColors[site,], line=-2.5, adj=0.1)

#   title(xlab = "Iteration")
#   title(ylab = "Parameter values (nM)")
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))), pos=ys[1], lwd=2,
#        labels=FALSE, tck=-0.01)
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
#   axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
#        labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=1, lwd=2, hadj=0.8)
#   axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
#        pos=1, lwd=2)
#   sapply(colnames(out.print), function(name) {
#           lines(x, out.print[, name], lwd=2, col=GetColorFunction(name))
#         }
#         )
#   lines(x, 10^probs.scaled, type="l", col="black")
#   dev.off()
# }

# PlotFlankKdRepression <- function(mirna,cutoff=FALSE,merge=FALSE,noncanon=TRUE){
#     par(kPlotParameters)

#   data<- GetSitesXCounts(mirna, "equilibrium", 5, 5, sitelist="paper")[,1,drop=FALSE]
#   kds <- GetKds(mirna,"equilibrium", 5, 5, "paper",nosite=TRUE)
#   data_seqs <- data[,1]
#   names(data_seqs) <- rownames(data)

#   print(data_seqs)

#   sites_all <- unlist(unique(subset(repression_df,mir==mirna,select=site_type)))
#   sites_all <- sites_all[sites_all!="nosite"]
#   site_seqs <- sapply(sites_all, function(site) {
#     data_seqs[site]
#   })
#   names(site_seqs) <- sites_all
#   print('sites_all')
#   print(sites_all)

#   print(site_seqs)
#   print(data_seqs)
#   site_seqs_noncanonical <- site_seqs
#   site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer"], site_seqs_noncanonical, invert=TRUE)]
#   site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer-A1"], site_seqs_noncanonical, invert=TRUE)]
#   site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer-m8"], site_seqs_noncanonical, invert=TRUE)]
#   site_seqs_noncanonical <- site_seqs_noncanonical[names(site_seqs_noncanonical)!="8mer-bG(6.7)"]
#   site_seqs_noncanonical <- site_seqs_noncanonical[names(site_seqs_noncanonical)!="7mer-m8bG(6.7)"]

#   site_seqs_canonical <- setdiff(sites_all,names(site_seqs_noncanonical))
#   print("canonical")
#   print(site_seqs_canonical)
#   print("noncanonical")
#   print(site_seqs_noncanonical)
#   if (merge == TRUE){
#     repression_df$site_type[which(repression_df$site_type %in% names(site_seqs_noncanonical))] <- "Noncanonical"
#       print(unique(repression_df$site_type))
#     sites_all <- c(site_seqs_canonical, "Noncanonical")
#   }

#   out <- sapply(sites_all, function(site) {
#     print(site)
#     if (cutoff != FALSE & site=="Noncanonical") {
#       reduced_frame <- subset(repression_df, mir==mirna & site_type==site & log_kd<=cutoff ,select=c(log_fc,log_kd))
#     } else {
#       reduced_frame <- subset(repression_df, 
#         mir==mirna & site_type==site,
#         select=c(log_fc,log_kd)
#       )      
#     }
#     print(dim(reduced_frame))
#     mean_values <- colMeans(reduced_frame)
#     sd_values <- apply(reduced_frame,2,sd)/sqrt(nrow(reduced_frame))
#     return(c(mean_values,sd_values))
#   })
#   par(mfrow=c(1,1))
#   colnames(out) <- sites_all
#   print(out)
#   if (noncanon == FALSE) {
#     out <- out[,c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8")]
#     sites_all <- colnames(out)
#   }


#   print(out)
#       xmin <- 1
#     xmax <- 2000
#     ymin <- -1
#     ymax <- 0.1
#     nosite_rep <- mean(subset(repression_df,mir==mirna & site_type=="nosite",select=c(log_fc))[,1])
#   plot(c(1,1),c(1,1),xlim=c(xmin, xmax), ylim = c(ymin,ymax),log='x',col="white",cex=2,pch=19, lwd=2,ann=FALSE,axes=FALSE)
#   arrows(1/(2^out[2,]), out[1,]+out[3,] - nosite_rep,1/(2^out[2,]),out[1,]-out[3,]-nosite_rep,length=0.05, lwd=1.5,angle=90, code=3)
#   arrows(1/(2^(out[2,]+out[4,])), out[1,]-nosite_rep,1/(2^(out[2,]-out[4,])),out[1,]-nosite_rep,length=0.05, lwd=1.5,angle=90, code=3)
#   points(1/(2^out[2,]),out[1,]-nosite_rep,col="black", bg=kSiteColors[sites_all,],pch=21,cex=1.5)


#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     ys <- seq(ymin,ymax,by=0.05)

#     xs <- xs[xs >= xmin & xs <= xmax]
#     ys <- ys[ys >= ymin & ys <= ymax]

#     # xmin <- min(xs)
#     # xmax <- max(xs)
#     # ymin <- min(ys)
#     # ymax <- max(ys)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#     yl <- seq(ymin,ymax, by= 0.1)





#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=ymin, lwd=0, cex.axis=1.7, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=ymin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=yl,
#          labels=round(yl,2),
#          pos=xmin, las=2, lwd=0, cex.axis=1.7,hadj=0.8)
#     axis(2, at=ys, labels=FALSE,
#          pos=xmin, lwd=2)


#   remove <- unique(which(is.na(out[1,])), which(is.na(out[1,])))
#   print(remove)
#   print(length(remove))
#   if (length(remove) > 0) {
#     out <- out[,-remove]  
#   }
#   text(x=800, y=0.05, mirna,cex=1.5)

#   text(x=800, y=-0.02, eval(substitute(expression(italic(r) == x), 
#             list(x = round(cor(out[2,],out[1,]),3)))),cex=1.5)
#   if (merge == TRUE) {
#   text(x=500, y=-0.09, eval(substitute(expression(log[2](italic(K)[D][italic(noncanon)]) <= x),
#             list(x = cutoff))), cex=1.5)
#   }

#   # Second plot
#   # kds <- GetKds(mirna,"equilibrium", 5, 5, "paper",nosite=TRUE)
#   # sites_all <- unlist(unique(subset(repression_df,mir==mirna,select=site_type)))
#   # sites_all <- sites_all[sites_all!="nosite"]
#   # out <- sapply(sites_all, function(site) {
#   #   reduced_frame <- subset(repression_df, mir==mirna & site_type==site,select=c(log_fc))
#   #   mean_values <- colMeans(reduced_frame)
#   #   sd_values <- apply(reduced_frame,2,sd)/sqrt(nrow(reduced_frame))
#   #   return(c(mean_values,log(kds[site],base=2),sd_values))
#   # })
#   # title(main=mirna,font.main=1,cex=2)
#   # plot(out[2,],out[1,],xlim=c(-10,1),ylim=c(-1,0.25),col=kSiteColors[sites_all,],cex=2,pch=seq(1,6),ann=FALSE,axes=FALSE)
#   # arrows(out[2,], out[1,]+out[3,],out[2,],out[1,]-out[3,],length=0.05, angle=90, code=3, lwd=1.5)
#   # points(out[2,],out[1,],col=kSiteColors[sites_all,],pch=seq(1,6),cex=2)
#   # axis(1,at=seq(-10,1),pos = -1,lwd = 2)
#   # axis(2, at = seq(-1,0.25,by=0.25),pos = -10,lwd=2)
#   out <<- out
#   sites_all <- sites_all[order(out[2,])]
#   legend(x=1.5,y=-0.5, legend=sites_all, bty="n", pch=19,col=kSiteColors[sites_all,],cex=1.1, ncol=1)
#   sites_all <<- sites_all
#   title(xlab=expression(italic(K)[D]))
#   title(ylab=expression(log[2](paste("fold change"))))
#   # text(x=30,y=0.1,round(cor(out[2,],out[1,]),3),cex=1.5)


# }

# PlotKdRepression <- function(mirna){
#     par(kPlotParameters)
#   kds <- GetKds(mirna,"equilibrium", 5, 5, "paper",nosite=TRUE)
#   sites_all <- unlist(unique(subset(repression_df,mir==mirna,select=site_type)))
#   sites_all <- sites_all[sites_all!="nosite"]
#   print(sites_all)
#   out <- sapply(sites_all, function(site) {
#     print(site)
#     reduced_frame <- subset(repression_df, mir==mirna & site_type==site,select=c(log_fc))
#     print(reduced_frame)
#     mean_values <- colMeans(reduced_frame)
#     sd_values <- apply(reduced_frame,2,sd)/sqrt(nrow(reduced_frame))
#     return(c(mean_values,log(kds[site],base=2),sd_values))
#   })
#   colnames(out) <- sites_all
#   print(out)
#   plot(out[2,],out[1,],xlim=c(-10,1),ylim=c(-1,0.25),col=kSiteColors[sites_all,],cex=2,pch=seq(1,5),ann=FALSE,axes=FALSE)
#   arrows(out[2,], out[1,]+out[3,],out[2,],out[1,]-out[3,],length=0.05, angle=90, code=3, lwd=1.5)
#   points(out[2,],out[1,],col=kSiteColors[sites_all,],pch=seq(1,5),cex=2)
#   axis(1,at=seq(-10,1),pos = -1,lwd = 2)
#   axis(2, at = seq(-1,0.25,by=0.25),pos = -10,lwd=2)
#   legend(x=, legend=sites_all,pch=seq(1,5),col=kSiteColors[sites_all,])

# }


# WriteIterationFile <- function(out,extension="") {
#   out.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                      "/equilibrium/kds/", site,"_flanking_",
#                      k.c.stockago, extension,".txt")
#   write.table(file=out.file, out, sep="\t", quote=FALSE, row.names=FALSE,
#                 col.names=TRUE)
# }

# WriteFinalParameterFile <- function(out,extension="") {
#   out.final <- out[dim(out)[1], ]
#   names(out.final) <- colnames(out)
#   out.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                            "/equilibrium/kds/final_", site,
#                            "_flanking_", k.c.stockago,extension, ".txt")
#   write.table(file=out.file, out.final, sep="\t", quote=FALSE, row.names=TRUE,
#               col.names=FALSE)
# }

# ## Helpful Fnctions

# CompareSiteInputsN <- function(mirna, experiment, n_constant_1, n_constant_2, sitelist,column) {
#   column_1 <- GetSitesXCounts(mirna, experiment, n_constant_1, sitelist)[,column,drop=FALSE]
#   column_2 <- GetSitesXCounts(mirna, experiment, n_constant_2, sitelist)[,column,drop=FALSE]
#   out <- cbind(column_1, column_2)
#   plot(out[,1], out[,2], log='xy')

#   return(out)
# }

# CompareSiteFlanksInputsN <- function(mirna, experiment, n_constant_1, n_constant_2, sitelist, site, column) {
#   column_1 <- GetSiteFlanksXCounts(mirna, experiment, n_constant_1, sitelist, site)[,column,drop=FALSE]
#   column_2 <- GetSiteFlanksXCounts(mirna, experiment, n_constant_2, sitelist, site)[,column,drop=FALSE]
#   out <- cbind(column_1, column_2)
#   plot(out[,1], out[,2], log='xy')
#   return(out)
# }










# MakeSiteIterationPlot <- function(out,method,mirna, num.kds, num.bgs, colors = FALSE) {
#   setEPS()
#   postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
#     mirna,"/", method, ".eps"))
#   par(kPlotParameters)
#   x = seq(1,dim(out)[1],length = max(ncol(out),1000))
#   out <- out[x,]
#   probs = out[ ,"-logp"]
#   out.print <- out[ , seq(dim(out)[2] - 1)]
#   out.print.kds <- Logistic(out[,1:num.kds],max = 1)
#   out.print.bgs <- 10^out.print[,(num.kds + 1) : (ncol(out.print))]
#   out.print <- out.print <- cbind(out.print.kds,out.print.bgs)
#   ys <- 10^c(max(floor(log10(min(out.print))), -5), ceiling(log10(max(out.print))))
#   probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])


#   plot(x    = x,
#        y    = 10^probs.scaled,
#        log  = "y",
#        axes = FALSE,
#        type = "l",
#        col = "white",
#        ylim = ys)
#   axis(side = 1,
#        at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
#        pos  = ys[1], labels = FALSE,
#        tck  = -0.01)
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
#   axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
#        labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=1, lwd=2)
#   title(main = mirna, font=1)
#   axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
#        pos=1)
#   if (colors != FALSE) {

#   cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
#   names(cols) <- colnames(out.print)
#   cols["AGO"] <- "grey"

#   sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=cols[name])
#     })
#   } else  {
#     sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=kSiteColors[name,])
#     })
#   }
#   lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
#   dev.off()
# }


# # Made FOR PAPER NOW
# GetSiteKds <- function(mirna, experiment, n_constant, sitelist) {
#     params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
#                  experiment, "/kds_PAPER/", n_constant, "_", 
#                  sitelist, "_PAPER.txt")
#     params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
#                          stringsAsFactors=FALSE))
#     return(params)
# }

# GetSiteKoffs <- function(mirna, experiment, n_constant, sitelist, costfunc, subset=FALSE, dil=FALSE) {
#     if (dil !=FALSE) {
#       dil.ext <- "_plusdil"
#     } else {
#       dil.ext <- ""
#     }
#     if (subset != FALSE) {
#       subset.ext <- paste0("_", subset, "-sites")
#     } else {
#       subset.ext <- ""
#     }
#     params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
#                  experiment, "/koffs_PAPER/", n_constant, "_", 
#                  sitelist, "_", costfunc, subset.ext, "_round2", dil.ext, ".txt")
#     print(params.file)
#     params <- data.frame(read.table(params.file, header=TRUE, row.names=NULL,
#                          stringsAsFactors=FALSE))
#     colnames(params) <- sapply(colnames(params),
#                                          RemovePrefixFromString, prefix="X")
#     params.final <- as.numeric(params[nrow(params), -ncol(params)])
#     names(params.final) <- colnames(params)[1:length(params.final)]
#     # print(params.final)
#     # kds <- params.final[grep("_Kd", names(params.final))]
#     # koffs <- params.final[grep("_koff", names(params.final))]
#     # kinetics_params <- cbind(kds, koffs)
#     # contam_params <- params.final["contam"]
#     return(params.final)
# }

# GetFlankKoffs <- function(mirna, experiment, n_constant, sitelist, site, costfunc, dil=FALSE) {
#     if (dil !=FALSE) {
#       dil.ext <- "_plusdil"
#     } else {
#       dil.ext <- ""
#     }
#     params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
#                  experiment, "/koffs_PAPER/", n_constant, "_", 
#                  sitelist, "_", costfunc, "_", site, "_round2", dil.ext, ".txt")
#     print(params.file)
#     params <- data.frame(read.table(params.file, header=TRUE, row.names=NULL,
#                          stringsAsFactors=FALSE))
#     colnames(params) <- sapply(colnames(params),
#                                          RemovePrefixFromString, prefix="X")
#     params.final <- as.numeric(params[nrow(params), -ncol(params)])
#     names(params.final) <- colnames(params)[1:length(params.final)]
#     # print(params.final)
#     # kds <- params.final[grep("_Kd", names(params.final))]
#     # koffs <- params.final[grep("_koff", names(params.final))]
#     # kinetics_params <- cbind(kds, koffs)
#     # contam_params <- params.final["contam"]
#     return(params.final)
# }

# CompareFlankingKoffs <- function(mirna, experiment, n_constant, sitelist, site1, site2, costfunc, dil=FALSE) {
#   flanks1 <- GetFlankKoffs(mirna, experiment, n_constant, sitelist, site1, costfunc, dil=FALSE)
#   flanks2 <- GetFlankKoffs(mirna, experiment, n_constant, sitelist, site2, costfunc, dil=FALSE)


# }

# PlotKoffsVsKds <- function(mirna, n_constant, sitelist, ident=FALSE, subset=FALSE, dil=FALSE, alpha=0.5) {
#     dev.new(xpos = 20, ypos = 20, height = 8, width = 12)
#     par(kPlotParameters)
#     par(mfrow = c(2, 3))
#     par(mgp = c(0.5, 0.6, 0.0))
#     costfunc <- "logres"
#     tempkds <- GetSiteKds(mirna, "equilibrium", n_constant, sitelist)
#     tempkoffs <- GetSiteKoffs(mirna, "kinetics", n_constant, sitelist, costfunc, subset=subset, dil = dil)
#     print(tempkoffs)
#     if (subset != FALSE) {
#       tempkds <- tempkds[c(1:subset, nrow(tempkds)-2), ]
#     } else {
#       tempkds <- tempkds[1:(nrow(tempkds)-2), ]
#     }
#     print(tempkoffs)
#       kinetics.kds <- 10^tempkoffs[1:nrow(tempkds)]    
#       kinetics.koff <- 10^tempkoffs[(nrow(tempkds) + 1):(2*nrow(tempkds))]
#       print(kinetics.kds)
#     tempkds <<- tempkds
#     tempkoffs <<- tempkoffs
#     kinetics.koff <<- kinetics.koff
#     kinetics.kds <<- kinetics.kds
#       break

#     if (sitelist %in% c("biological12mers", "biological12mersnew")){
#       colors <- rep("black", length(kinetics.kds))
#       ind.6mers <- grep(GiveTargetSequence("miR-1", 2, 7), rownames(tempkds), fixed=TRUE)
#       colors[ind.6mers] <- ConvertRColortoRGB(kSiteColors["6mer", ], alpha = alpha)
#       ind.7mersA1 <- grep(GiveTargetSequence("miR-1", 1, 7), rownames(tempkds), fixed=TRUE)
#       colors[ind.7mersA1] <- ConvertRColortoRGB(kSiteColors["7mer-A1", ], alpha = alpha)
#       ind.7mersm8 <- grep(GiveTargetSequence("miR-1", 2, 8), rownames(tempkds), fixed=TRUE)
#       colors[ind.7mersm8] <- ConvertRColortoRGB(kSiteColors["7mer-m8", ], alpha = alpha)
#       ind.8mers <- grep(GiveTargetSequence("miR-1", 1, 8), rownames(tempkds), fixed=TRUE)
#       colors[ind.8mers] <- ConvertRColortoRGB(kSiteColors["8mer", ], alpha = alpha)
#     } else {
#       colors <- kSiteColors[rownames(tempkds), ]
#     }

#     kdl    <- seq(-4, 3)
#     kds    <- sapply(kdl, function(i) {10^i * seq(9)})
#     kdl <- 10^kdl
#     kdmin <- kdl[1]
#     kdmax <- kdl[length(kdl)]
#     kdlims <- c(kdmin, kdmax)
#     koffl    <- seq(-3, 2)
#     koffs    <- sapply(koffl, function(i) {10^i * seq(9)})
#     koffl    <- 10^koffl
#     koffmin <- koffl[1]
#     koffmax <- koffl[length(koffl)]
#     kofflims <- c(koffmin, koffmax)

#     plot(tempkds$Mean,
#       kinetics.kds,
#       log = 'xy',
#       col=colors,
#       xlim = kdlims,
#       ylim = kdlims,
#       ann=FALSE,
#       axes=FALSE)
#     mtext(costfunc, side = 3, line = 1)
#     abline(0, 1, lty = 2, lwd = 0.5)
#     # axis(1, at = ticks, labels = FALSE, pos = min(xlims))
#     # axis(2, at = ticks, labels = FALSE, pos = min(xlims))
#     # axis(1, at = 10^seq(-3, 3), lwd = 0, pos = min(xlims))
#     # axis(2, at = 10^seq(-3, 3), lwd = 0, pos = min(xlims))
#     axis(1, at=kdl,
#          labels=sapply(kdl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=kdmin, lwd=0, line = 3, cex.axis = 1.5)
#     axis(1, at=kds, labels=FALSE,
#          pos=kdmin)
#     # Label the axis at each order of magnitude.

#     axis(2, at=kdl,
#          labels=sapply(kdl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
#     axis(2, at=kds, labels=FALSE,
#          pos=kdmin)


#     title(xlab = "Equilibrium KD (nM)", line= 1, cex.lab = 1.5)
#     title(ylab = "Kinetics KD (nM)", line = 1.5, cex.lab = 1.5)



#     plot(tempkds$Mean,
#       kinetics.koff,
#       log = 'xy',
#       col=colors,
#       xlim = kdlims,
#       ylim = kofflims,
#       ann = FALSE,
#       axes = FALSE)
#     xstart <- kdmin
#     range <- kdmax/kdmin
#     yoffset <- (kinetics.koff["None_koff"])/(tempkds$Mean[nrow(tempkds)])
#     segments(xstart, xstart*yoffset, xstart*range, xstart*yoffset*range, lty = 2, lwd = 0.5)

#     axis(1, at=kdl,
#          labels=sapply(kdl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=koffmin, lwd=0, line = 3, cex.axis = 1.5)
#     axis(1, at=kds, labels=FALSE,
#          pos=koffmin)
#     # Label the axis at each order of magnitude.

#     axis(2, at=koffl,
#          labels=sapply(koffl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
#     axis(2, at=koffs, labels=FALSE,
#          pos=kdmin)


#     title(xlab = "Equilibrium KD (nM)", line= 1, cex.lab = 1.5)
#     title(ylab = "Kinetics koff (min^-1)", line = 1.5, cex.lab = 1.5)





#     if (ident == TRUE) {
#       identify(tempkds$Mean[1:(nrow(tempkds)-2)],
#       kinetics.koff,
#       labels = rownames(tempkds)[1:(nrow(tempkds)-2)])
#     }
#     title(main = mirna)
#     plot(kinetics.kds,
#       kinetics.koff,
#       log = 'xy',
#       col=colors,
#       xlim = kdlims,
#       ylim = kofflims,
#       ann = FALSE,
#       axes = FALSE)
#     xstart <- kdmin
#     range <- kdmax/kdmin
#     print(kinetics.koff)
#     yoffset <- (kinetics.koff["None_koff"])/(kinetics.kds["None_Kd"])
#     segments(xstart, xstart*yoffset, xstart*range, xstart*yoffset*range, lty = 2, lwd = 0.5)

#     axis(1, at=kdl,
#          labels=sapply(kdl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=koffmin, lwd=0, line = 3, cex.axis = 1.5)
#     axis(1, at=kds, labels=FALSE,
#          pos=koffmin)
#     # Label the axis at each order of magnitude.

#     axis(2, at=koffl,
#          labels=sapply(koffl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
#     axis(2, at=koffs, labels=FALSE,
#          pos=kdmin)


#     title(xlab = "Kinetics KD (nM)", line= 1, cex.lab = 1.5)
#     title(ylab = "Kinetics koff (min^-1)", line = 1.5, cex.lab = 1.5)






#     if (ident == TRUE) {
#       identify(tempkds$Mean[1:(nrow(tempkds)-2)],
#       kinetics.koff,
#       labels = rownames(tempkds)[1:(nrow(tempkds)-2)])
#     }
#     costfunc <- "multinom"
#     tempkoffs <- GetSiteKoffs(mirna, "kinetics", n_constant, sitelist, costfunc, subset=subset, dil = dil)
#     print(tempkoffs)
#       kinetics.kds <- 10^tempkoffs[1:nrow(tempkds)]
#       print(kinetics.kds)
#       break
#       kinetics.koff <- 10^tempkoffs[(nrow(tempkds) + 1):(2*nrow(tempkds))]
#     tempkds <<- tempkds
#     tempkoffs <<- tempkoffs
#     kinetics.koff <<- kinetics.koff
#     kinetics.kds <<- kinetics.kds
#     if (sitelist %in% c("biological12mers", "biological12mersnew")){
#       colors <- rep("black", length(kinetics.kds))
#       ind.6mers <- grep(GiveTargetSequence("miR-1", 2, 7), rownames(tempkds), fixed=TRUE)
#       colors[ind.6mers] <- ConvertRColortoRGB(kSiteColors["6mer", ], alpha = alpha)
#       ind.7mersA1 <- grep(GiveTargetSequence("miR-1", 1, 7), rownames(tempkds), fixed=TRUE)
#       colors[ind.7mersA1] <- ConvertRColortoRGB(kSiteColors["7mer-A1", ], alpha = alpha)
#       ind.7mersm8 <- grep(GiveTargetSequence("miR-1", 2, 8), rownames(tempkds), fixed=TRUE)
#       colors[ind.7mersm8] <- ConvertRColortoRGB(kSiteColors["7mer-m8", ], alpha = alpha)
#       ind.8mers <- grep(GiveTargetSequence("miR-1", 1, 8), rownames(tempkds), fixed=TRUE)
#       colors[ind.8mers] <- ConvertRColortoRGB(kSiteColors["8mer", ], alpha = alpha)
#     } else {
#       colors <- kSiteColors[rownames(tempkds), ]
#     }

#     plot(tempkds$Mean,
#       kinetics.kds,
#       log = 'xy',
#       col=colors,
#       xlim = kdlims,
#       ylim = kdlims,
#       ann=FALSE,
#       axes=FALSE)
#     mtext(costfunc, side = 3, line = 1)
#     abline(0, 1, lty = 2, lwd = 0.5)
#     # axis(1, at = ticks, labels = FALSE, pos = min(xlims))
#     # axis(2, at = ticks, labels = FALSE, pos = min(xlims))
#     # axis(1, at = 10^seq(-3, 3), lwd = 0, pos = min(xlims))
#     # axis(2, at = 10^seq(-3, 3), lwd = 0, pos = min(xlims))
#     axis(1, at=kdl,
#          labels=sapply(kdl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=kdmin, lwd=0, line = 3, cex.axis = 1.5)
#     axis(1, at=kds, labels=FALSE,
#          pos=kdmin)
#     # Label the axis at each order of magnitude.

#     axis(2, at=kdl,
#          labels=sapply(kdl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
#     axis(2, at=kds, labels=FALSE,
#          pos=kdmin)


#     title(xlab = "Equilibrium KD (nM)", line= 1, cex.lab = 1.5)
#     title(ylab = "Kinetics KD (nM)", line = 1.5, cex.lab = 1.5)



#     plot(tempkds$Mean,
#       kinetics.koff,
#       log = 'xy',
#       col=colors,
#       xlim = kdlims,
#       ylim = kofflims,
#       ann = FALSE,
#       axes = FALSE)
#     xstart <- kdmin
#     range <- kdmax/kdmin
#     yoffset <- (kinetics.koff["None_koff"])/(tempkds$Mean[nrow(tempkds)])
#     segments(xstart, xstart*yoffset, xstart*range, xstart*yoffset*range, lty = 2, lwd = 0.5)

#     axis(1, at=kdl,
#          labels=sapply(kdl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=koffmin, lwd=0, line = 3, cex.axis = 1.5)
#     axis(1, at=kds, labels=FALSE,
#          pos=koffmin)
#     # Label the axis at each order of magnitude.

#     axis(2, at=koffl,
#          labels=sapply(koffl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
#     axis(2, at=koffs, labels=FALSE,
#          pos=kdmin)


#     title(xlab = "Equilibrium KD (nM)", line= 1, cex.lab = 1.5)
#     title(ylab = "Kinetics koff (min^-1)", line = 1.5, cex.lab = 1.5)





#     if (ident == TRUE) {
#       identify(tempkds$Mean[1:(nrow(tempkds)-2)],
#       kinetics.koff,
#       labels = rownames(tempkds)[1:(nrow(tempkds)-2)])
#     }
#     title(main = mirna)
#     plot(kinetics.kds,
#       kinetics.koff,
#       log = 'xy',
#       col=colors,
#       xlim = kdlims,
#       ylim = kofflims,
#       ann = FALSE,
#       axes = FALSE)
#     xstart <- kdmin
#     range <- kdmax/kdmin
#     print(kinetics.koff)
#     yoffset <- (kinetics.koff["None_koff"])/(kinetics.kds["None_Kd"])
#     segments(xstart, xstart*yoffset, xstart*range, xstart*yoffset*range, lty = 2, lwd = 0.5)

#     axis(1, at=kdl,
#          labels=sapply(kdl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=koffmin, lwd=0, line = 3, cex.axis = 1.5)
#     axis(1, at=kds, labels=FALSE,
#          pos=koffmin)
#     # Label the axis at each order of magnitude.

#     axis(2, at=koffl,
#          labels=sapply(koffl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
#     axis(2, at=koffs, labels=FALSE,
#          pos=kdmin)


#     title(xlab = "Kinetics KD (nM)", line= 1, cex.lab = 1.5)
#     title(ylab = "Kinetics koff (min^-1)", line = 1.5, cex.lab = 1.5)






#     if (ident == TRUE) {
#       identify(tempkds$Mean[1:(nrow(tempkds)-2)],
#       kinetics.koff,
#       labels = rownames(tempkds)[1:(nrow(tempkds)-2)])
#     }

# }


# # PlotKoffsVsKds("miR-1", 5, "biological12mersnew", "logres", dil=TRUE, alpha=0.3)
# # break

# GetFlankbgKds <- function(mirna, experiment, n_constant, sitelist, site) {
#     params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
#                  experiment, "/kds_PAPER/", n_constant, "_", 
#                  sitelist, "_", site, "_controlplfold_PAPER.txt")

#     params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
#                          stringsAsFactors=FALSE))
#     return(params)
# }


# GetFlankbyPlFoldKds <- function(mirna, experiment, n_constant, sitelist, site) {
#     params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
#                  experiment, "/kds_PAPER/", n_constant, "_", 
#                  sitelist, "_", site, "_withprob_PAPER.txt")

#     params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
#                          stringsAsFactors=FALSE))
#     return(params)
# }

# PlotFlankPlFoldKds <- function(mirna, experiment, n_constant, sitelist, site) {
#   kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
#   kds.plot <- kds$Mean
#   kds.sub <- GetFlankbyPlFoldKds(mirna, experiment, n_constant, sitelist, site)
#   sub.divisions <- unique(sapply(rownames(kds.sub), function(name) {
#     temp_name <- paste0("_",unlist(strsplit(name, split = "_"))[2], collapse = "")
#     }))
#   print(sub.divisions)
#   par(kPlotParameters)
#   average_vector <- matrix(NaN,nrow=length(sub.divisions), ncol=length(kds.plot))
#   plot(kds.plot, rnorm(length(kds.plot),mean=0,sd=0.1),
#        log='x', xlim = c(0.0001, 10), ylim = c(0, length(sub.divisions)+5),
#        col=sapply(rownames(kds), GetColorFunction))
#   starting_sd <- sd(log10(kds.plot))^2
#   text(10,0,sd(log10(kds.plot))^2/starting_sd)
#   sapply(1:length(sub.divisions), function(ind) {
#     print(ind)
#     print(sub.divisions)
#     print(sub.divisions[ind])
#     print(paste0(sub.divisions[ind],"$",collapse=""))
#     sub.division <- paste0(sub.divisions[ind],"$",collapse="")
#     print(sub.division)
#     kds.plot <- kds.sub[grep(sub.division,rownames(kds.sub),perl=TRUE),]$Median
#     print(length(kds.plot))
#     print(length(log10(kds.plot)-mean(log10(kds.plot))))
#     print(length(average_vector[ind,]))
#     average_vector[ind,] <<- log10(kds.plot)-mean(log10(kds.plot))+mean(log10(kds$Mean))
#     kd.plot_alt <- 10^(log10(kds.plot)-mean(log10(kds.plot))+mean(log10(kds$Mean)))
#     points(kds.plot, ind+rnorm(length(kds.plot),mean=0,sd=0.1), col=sapply(rownames(kds), GetColorFunction))
#      text(10,ind,round(sd(log10(kds.plot))^2/starting_sd,2))
#   })
#   final <- 10^colMeans(average_vector)
#   points(10^colMeans(average_vector),
#          length(sub.divisions)+2+rnorm(length(kds.plot), mean=0, sd=0.1),
#          col=sapply(rownames(kds), GetColorFunction))
#        text(10,length(sub.divisions)+2,round(sd(log10(final))^2/starting_sd,2))

# }

# PlotFlankControlKds <- function(mirna, experiment, n_constant, sitelist, site) {
#   kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
#   kds.plot <- kds$Mean
#   kds.sub <- GetFlankbgKds(mirna, experiment, n_constant, sitelist, site)$Mean

#   par(kPlotParameters)
#   # par(mfrow=c(2,1))
#   # plot(kds.plot, 1+rnorm(length(kds.plot),mean=0,sd=0.1),
#   #      log='x', xlim = c(0.0001, 10), ylim = c(0, 5),
#   #      col=sapply(rownames(kds), GetColorFunction))

#   # points(kds.sub, 2 + rnorm(length(kds.plot),mean=0,sd=0.1),
#   #        col=sapply(rownames(kds), GetColorFunction))
#       xmin = 0.0001
#       xmax = 10
#         plot( kds.sub, kds.plot,
#        log='xy', xlim = c(0.0001, 10), ylim = c(0.0001, 10),
#        col=sapply(rownames(kds), GetColorFunction))
#       linear_model <- lm(log10(kds.plot) ~ log10(kds.sub))
#       m <- linear_model$coefficients[2]
#       b <- linear_model$coefficients[1]

#     x_line <- 10^seq(log10(xmin), log10(xmax),length = 20)
#     y_line <- 10^(m*log10(x_line) + b)

#     lines(x_line, y_line, lty = 2,lwd = 0.5)
#     text(1e-1, 1e-3, bquote(log(y) == log(x)*.(round(m,2)) + .(round(b,2))))
#     text(1e-1, 2e-3, round(cor(log(kds.sub), log(kds.plot))^2,2))
#     # text(x = 0.001, y = 0.95, bquote(italic(r)^2 ==  .(cor_text)))


# }


# PlotSiteKds <- function(mirna, experiment, n_constant, sitelist, plotlist,
#                         xpos = 20, ypos = 20) {
#   # Get kds for all site-types of the mirna.
#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)

#   kds <- GetSiteKds(mirna, experiment, n_constant, sitelist)
#   if (length(plotlist) == 0) {
#     site_list <- rownames(data)
#   } else if (class(plotlist) == "character") {
#     site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
#                                             "computation/AgoRBNS/",
#                                             "AssignSiteTypes/sites.", mirna,
#                                             "_", plotlist, ".txt"),
#                               stringsAsFactors=FALSE)[,1], "None")
#   } else {
#     site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
#   }
#   kds <- kds[-which(rownames(kds) %in% c("bg", "AGO")),]
#   kds <- kds[order(kds$Full),]

#   par(kPlotParameters)
#   xs <- kds
#   ys <- nrow(kds) - seq(nrow(kds)) + 1
#   plot(1, type ="n",
#           axes    = FALSE,
#           log = 'x',
#           ylim       = c(0, 22),
#           xlim       = rev(c(0.00003, 3)))
#  arrows(kds$Upper_CI, nrow(kds) - seq(nrow(kds)) + 1,
#         kds$Lower_CI, nrow(kds) - seq(nrow(kds)) + 1, length=0.05, angle=90, code=3)

#   title(main = mirna,
#         line = -2,
#         adj  = 0.1)
#   title(xlab = expression(italic(K)[D]))
#  points(kds$Mean,nrow(kds) - seq(nrow(kds)) + 1,
#           col = "black",
#           bg = kSiteCategoryColors[rownames(kds),],
#           pch = 21,
#           lwd=1,
#           cex = 1.2)
#   ymin=0.0001
#   ymax=3
#   ys <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
#   ys <- ys[which(ys>=ymin & ys <= ymax)]
#   # ys <- ys[ys >= ymin & ys <= ymax]
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#   axis(1, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=0, lwd=0)
#   axis(1, at=ys, labels=FALSE,
#        pos=0, lwd=1)
# # points(rep(0.0006, nrow(kds)),nrow(kds) - seq(nrow(kds)) + 1,col=kSiteCategoryColors[rownames(kds),],cex=1.5,pch=19)


#   text((kds$Full)/1.6,nrow(kds) - seq(nrow(kds)) + 1, labels=ConvertTtoU(rownames(kds)), adj=0, cex = 0.8, col= "black")


#   if (mirna == "miR-124") {
#       legend(x = 10^-2.1, y = 8,legend = c("7-8-nt canonical site",
#                                  "6-nt canonical site",
#                                  "Enhanced 6mer-containing sites",
#                                  "Noncanonical sites",
#                                  "3' sites"),
#          bty="n",
#          col="black",
#          pt.bg = c("purple2", "deepskyblue2", "cyan", "violet", "green3"),
#          pch = 21,
#          pt.cex = 1.2,
#          pt.lwd=1)

#       } else if (mirna == "miR-155") {
#   legend(x = 10^-2.1, y = 8,legend = c("7-8-nt canonical site",
#                                  "6-nt canonical site",
#                                  "Noncanonical sites",
#                                  "3' sites",
#                                  "??"),
#          bty="n",
#          col="black",
#          pt.bg = c("purple2", "deepskyblue2", "violet", "green3","gray"),
#          pch = 21,
#          pt.cex = 1.2,
#          pt.lwd=1)
# } else if (mirna == "miR-1") {
#     legend(x = 10^-2.1, y = 7,legend = c("7-8-nt canonical site",
#                                  "6-nt canonical site",
#                                  "Noncanonical sites",
#                                  "??"),
#          bty="n",
#          col="black",
#          pt.bg = c("purple2", "deepskyblue2", "violet", "gray"),
#          pch = 21,
#          pt.cex = 1.2,
#          pt.lwd=1)

# } else if (mirna == "let-7a") {
#     legend(x = 10^-2.1, y = 6.5,legend = c("7-8-nt canonical site",
#                                  "6-nt canonical site",
#                                  "Noncanonical sites",
#                                  "??"),
#          bty="n",
#          col="black",
#          pt.bg = c("purple2", "deepskyblue2", "violet", "gray"),
#          pch = 21,
#          pt.cex = 1.2,
#          pt.lwd=1)

# } else {
#   legend(x = 10^-2.1, y = 9.5,legend = c("7-8-nt canonical site",
#                                  "6-nt canonical site",
#                                  "Enhanced 6mer-containing sites",
#                                  "Noncanonical sites",
#                                  "3' sites",
#                                  "??"),
#          bty="n",
#          col="black",
#          pt.bg = c("purple2", "deepskyblue2", "cyan", "violet", "green3","gray"),
#          pch = 21,
#          pt.cex = 1.2,
#          pt.lwd=1)

#   }
# }

# PlotSiteEnrichments <- function(mirna, experiment, n_constant, sitelist,
#                                 plotlist, xpos = 20, ypos = 20,
#                                 bgoff = FALSE) {
#   params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
#   print(params)
#   sitesXcounts <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
#   sitesXcounts <- sitesXcounts[,-1]
#   print(sitesXcounts)
#   # Remove any rows with no reads whoatsoever.
#   sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
#   # Remove any rows for which there are no input reads:
#   sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
#   kds <- params[1:(nrow(params)-2),]$Mean
#   names(kds) <- rownames(sitesXcounts)

#   bgs <- rep(params["bg",]$Mean, 5)
#   k.c.stockago <- params["AGO",]$Mean
#   k.c.stockago.plot <- stockago[mirna, "equilibrium"]
#   c.I.tots <- Norm(sitesXcounts[,2])*100
#   names(c.I.tots) <- rownames(sitesXcounts)
#   data <- sitesXcounts[,3:7]
#   names(c.I.tots) <- rownames(data)
#   x_model <- seq(min(as.numeric(colnames(data))) * 0.7,
#                  max(as.numeric(colnames(data))) / 0.7,
#                  length=100)
#   x_model.alt <- colnames(data)
#   c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
#   c.totals.alt <<- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)

#   colnames(c.totals) <- colnames(x_model)
#   rownames(c.totals) <- rownames(data)

#   c.agos <- sapply(x_model, function(x) {
#       as.numeric(x) * k.c.stockago / 100
#     }
#   )
#   print(k.c.stockago)
#   c.agos.alt <<- sapply(x_model.alt, function(x) {
#       as.numeric(x) * k.c.stockago / 100
#     }
#   )

#   c.bounds <- as.matrix(
#     sapply(c.agos, function(x) {
#         return(GetBoundRNA(kds, c.I.tots, x))
#       }
#     )
#   )

#     c.bounds.alt <<- as.matrix(
#     sapply(c.agos.alt, function(x) {
#       print(x)
#         return(GetBoundRNA(kds, c.I.tots, x))
#       }
#     )
#   )

#   c.frees <- c.totals - c.bounds
#   c.frees.alt <<- c.totals.alt - c.bounds.alt
#   c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
#   c.bgs.alt <<- t(t(c.frees.alt) * bgs / colSums(c.frees.alt))

#   if (bgoff == TRUE) {
#     c.all <- c.bounds
#     c.all.alt <- c.bounds.alt

#   } else {
#     c.all <- c.bounds + c.bgs
#     c.all.alt <- c.bounds.alt + c.bgs.alt

#   }

#   c.final <- data.frame(t(t(c.all) / colSums(c.all)))
#   c.final.alt <<- data.frame(t(t(c.all.alt) / colSums(c.all.alt)))

#   rownames(c.final) <- rownames(data)
#   rownames(c.final.alt) <- rownames(data)

#   x <- c(40,12.65,4,1.265,0.4)*k.c.stockago.plot/100*1000
#   y <- c(1,1,1,1,1)
#   sites.norm <<- Norm(c.I.tots)
#   data.norm <- t(t(data)/colSums(data))

#   data.R <<- data.norm/(sites.norm)
#   model.R <<- c.final/(sites.norm)
#   model.R.alt <<- c.final.alt/(sites.norm)

#   xmin <- min(x)*0.3
#   xmax <- max(x)*3
#   ymin <- 0.2
#   ymax <- 300
#   yextension <- (ymax/ymin)

#   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   xs <- xs[xs >= xmin & xs <= xmax]
#   xmin <- min(xs)
#   xmax <- max(xs)
#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))
#   ys <- ys[ys >= ymin & ys <= ymax]
#   ymin <- min(ys)
#   ymax <- max(ys)
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 6)
#   par(kPlotParameters)

#   plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.4), ylim=c(ymin, ymax), type="l",
#      col="white", axes=FALSE, ann=FALSE)        
#   # Generate tickmarks for axis.

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=xmin, las=2, lwd=0)
#   axis(2, at=ys, labels=FALSE,
#        pos=xmin)

#   title(main = mirna, font.main=1, line=-2, adj=0.1)
#   title(xlab = "[AGO2-miRNA] (pM)", adj=0.3)
#   title(ylab = "Enrichment")
#   mirna.trim = paste0(strsplit(mirna, split = "-")[[1]][1:2],collapse="-")
#   print(mirna.trim)
#   plotlist <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
#                                           "computation/AgoRBNS/",
#                                           "AssignSiteTypes/sites.", mirna.trim,
#                                             "_", plotlist, ".txt"),
#                               stringsAsFactors=FALSE)[,1], "None")

#   legend.names <- rownames(data)[order(kds)]
#   legend.names <- legend.names[which(legend.names %in% plotlist)]
#   ordered_list <- legend.names


#   legend(x=xmax, y=ymax, legend=ConvertTtoU(legend.names), pch=19,
#          col=kSiteColors[ordered_list], bty="n")
#   for (name in plotlist) {
#     type = "p"
#     points(x, data.R[name, ], col=kSiteColors[name], type=type, pch=19, cex=1.2)
#     lines(x_model*k.c.stockago.plot/100*1000, model.R[name, ], col=kSiteColors[name],lwd=2)      
#   }
# }

# GetEquilibriumModel <- function(mirna, experiment, n_constant, sitelist) {
#   # Get the relevant parameters:
#   params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
#   # Get the experimental sequencing counts:
#   sXc <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
#   # Remove the "Seq" column
#   sXc <- sXc[,-1]

#   # Remove any rows with no reads whoatsoever.
#   sXc <- sXc[rowSums(sXc)>0,]
#   # Remove any rows for which there are no input reads:
#   sXc <- sXc[sXc[,1]>0,]
#   kds <- params[1:(nrow(params)-2),]$Mean
#   names(kds) <- rownames(sXc)

#   bgs <- rep(params["bg",]$Mean, 5)
#   k.c.stockago <- params["AGO",]$Mean

#   c.I.tots <- Norm(sXc[,2])*100
#   names(c.I.tots) <- rownames(sXc)
#   data <- sXc[,3:7]
#   names(c.I.tots) <- rownames(data)
#   x_points <- c(40,12.65,4,1.265,0.4)

#   c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = length(x_points), byrow=FALSE)
#   colnames(c.totals) <- colnames(x_points)
#   rownames(c.totals) <- rownames(data)

#   c.agos <- sapply(x_points, function(x) {
#       as.numeric(x) * k.c.stockago / 100
#     }
#   )
#   names(c.agos) <- x_points

#   c.bounds <- as.matrix(
#     sapply(c.agos, function(x) {
#         return(GetBoundRNA(kds, c.I.tots, x))
#       }
#     )
#   )
#   return(c.bounds / c.totals)

# }



# PlotSiteOccupancy <- function(mirna, experiment, n_constant, sitelist,
#                                 plotlist, xpos = 20, ypos = 20,
#                                 bgoff = FALSE) {
#   params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
  
#   sitesXcounts <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
#   sitesXcounts <- sitesXcounts[,-1]

#   # Remove any rows with no reads whoatsoever.
#   sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
#   # Remove any rows for which there are no input reads:
#   sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
#   kds <- params[1:(nrow(params)-2),]$Mean
#   names(kds) <- rownames(sitesXcounts)

#   bgs <- rep(params["bg",]$Mean, 5)
#   k.c.stockago <- params["AGO",]$Mean
#   k.c.stockago.plot <- stockago[mirna, "equilibrium"]

#   c.I.tots <- Norm(sitesXcounts[,2])*100
#   names(c.I.tots) <- rownames(sitesXcounts)
#   data <- sitesXcounts[,3:7]
#   names(c.I.tots) <- rownames(data)
#   x_model <- seq(min(as.numeric(colnames(data))) * 0.7,
#                  max(as.numeric(colnames(data))) / 0.7,
#                  length=100)
#   x_points <- c(40,12.65,4,1.265,0.4)

#   c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
#   colnames(c.totals) <- colnames(x_model)
#   rownames(c.totals) <- rownames(data)

#   c.agos_model <- sapply(x_model, function(x) {
#       as.numeric(x) * k.c.stockago / 100
#     }
#   )

#   c.agos_points <- sapply(x_points, function(x) {
#       as.numeric(x) * k.c.stockago / 100
#     }
#   )
#   c.bounds_model <- as.matrix(
#     sapply(c.agos_model, function(x) {
#         return(GetBoundRNA(kds, c.I.tots, x))
#       }
#     )
#   )

#   c.bounds_points <- as.matrix(
#     sapply(c.agos_points, function(x) {
#         return(GetBoundRNA(kds, c.I.tots, x))
#       }
#     )
#   )

#   c.bounds_model_norm <- apply(c.bounds_model, 2, Norm)
#   c.bounds_points_norm <- apply(c.bounds_points, 2, Norm)
#   rownames(c.bounds_model_norm) <- rownames(sitesXcounts)
#   rownames(c.bounds_points_norm) <- rownames(sitesXcounts)
#   xmin <- min(x_points*k.c.stockago.plot/100*1000)*0.3
#   xmax <- max(x_points*k.c.stockago.plot/100*1000)*3
#   ymin <- 0
#   ymax <- 0.7
#   yextension <- (ymax/ymin)
#   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   xs <- xs[xs >= xmin & xs <= xmax]
#   xmin <- min(xs)
#   xmax <- max(xs)
#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   ys <- seq(0,7)/10
#   yl <- ys
#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 6)
#   par(kPlotParameters)

#   plot(1,
#        type = "n",
#        log='x',
#        xlim=c(xmin, xmax*(xmax/xmin)^0.4),
#        ylim=c(ymin, ymax),
#        axes=FALSE,
#        ann=FALSE)        
#   # Generate tickmarks for axis.

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=yl,
#        pos=xmin, las=2, lwd=0)
#   axis(2, at=ys, labels=FALSE,
#        pos=xmin)

#   title(main = mirna, font.main=1, line=-2, adj=0.1)
#   title(xlab = "[AGO2-miRNA] (pM)", adj=0.3)
#   title(ylab = "Enrichment")

#   plotlist <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
#                                           "computation/AgoRBNS/",
#                                           "AssignSiteTypes/sites.", mirna,
#                                             "_", plotlist, ".txt"),
#                               stringsAsFactors=FALSE)[,1], "None")

#   legend.names <- rownames(data)[order(kds)]
#   legend.names <- legend.names[which(legend.names %in% plotlist)]
#   ordered_list <- legend.names


#   legend(x=xmax, y=ymax, legend=ConvertTtoU(legend.names), pch=19,
#          col=kSiteColors[ordered_list], bty="n")
#   for (name in plotlist) {
#     type = "p"
#     points(x_points*k.c.stockago.plot/100*1000, c.bounds_points_norm[name,], col=kSiteColors[name], type=type, pch=19,cex=1.2)
#     lines(x_model*k.c.stockago.plot/100*1000, c.bounds_model_norm[name,], col=kSiteColors[name], lwd = 2)      
#   }
# }



# PlotSiteScatterWithInput <- function(mirna, experiment, n_constant, sitelist, column,
#                                      xpos = 20, ypos = 20) {
  
#   sitesXcounts <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
#   sitesXcounts <- sitesXcounts[,-1]

#   # Remove any rows with no reads whoatsoever.
#   sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
#   # Remove any rows for which there are no input reads:
#   sitesXcounts <- sitesXcounts[sitesXcounts[,2]>0,]

#   x <- Norm(sitesXcounts[,2])*100
#   y <- Norm(sitesXcounts[,column])*100

#   xymin <- 0.05
#   xymax <- 100
#   yextension <- (xymax/xymin)
#   xys <- c(sapply(seq(floor(log10(xymin)), ceiling(log10(xymax))), function(x) seq(10)*10^x))
#   xys <- xys[xys >= xymin & xys <= xymax]
#   xymin <- min(xys)
#   xymax <- max(xys)
#   xyl <- 10^seq(ceiling(log10(min(xys))), floor(log10(max(xys))))

#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 6)
#   par(kPlotParameters)
#   # Make dummy plot:
#   plot(1, type = "n",log='xy',
#        xlim=c(xymin, xymax*(xymax/xymin)^0.4),
#        ylim=c(xymin, xymax),
#       axes=FALSE, ann=FALSE)
#   # Make x=y line:
#   segments(xymin, xymin, xymax, xymax, lty = 2)
#   # Make the lines connecting the points to the x = y line:
#   segments(x,x,x,y,lt = 2, col = kSiteColors[rownames(sitesXcounts), ])
#   # Make axes:
#   axis(1, at=xyl,
#        labels=xyl,
#        pos=xymin, lwd=0)
#   axis(1, at=xys, labels=FALSE,
#        pos=xymin)
#   axis(2, at=xyl,
#        labels=xyl,
#        pos=xymin, las=2, lwd=0)
#   axis(2, at=xys, labels=FALSE,
#        pos=xymin)
#   # Add the points to the plot:
#   points(x, y,
#      col=kSiteColors[rownames(sitesXcounts),])        

#   ago.percent <- as.numeric(colnames(sitesXcounts)[column])/100
#   title(main = paste0(ago.percent*stockago[mirna,"equilibrium"]*1000,' pM AGO2-', mirna), line=-2.5, adj=0.1)
#   title(xlab = "Input library (%)", adj=0.4)
#   title(ylab = "AGO-bound library(%)")


#   legend(x=50, y=2, legend=ConvertTtoU(rownames(sitesXcounts)), pch=19,
#          col=kSiteColors[rownames(sitesXcounts), ], bty="n", y.intersp=0.9)
#   }




# PlotSiteFlankEnrichments <- function(mirna, experiment, n_constant, sitelist, plotlist, site, bgoff = FALSE, xpos = 20, ypos = 20) {
#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 6)
#   par(kPlotParameters)
#   params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
#   flank.kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
#   sitesXcounts <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
#   sitesXcounts <- sitesXcounts[,-1]

#   s.c <- as.numeric(sitesXcounts[site, ])

#   sfXc <- GetSiteFlanksXCounts(mirna, experiment, n_constant, sitelist, site)

#   colnames(sfXc)[1] <- "I"
#   sfXc <- sfXc[rowSums(sfXc[, 2:6]) > 0,]
#   sfXc <- sfXc[which(sfXc[,1] > 0),,drop=FALSE]

#   sfXc <- round(t((t(sfXc) / colSums(sfXc)) * s.c))
#   sfXc[is.na(sfXc)] <- 0

#   # Remove any rows with no reads whoatsoever.
#   sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
#   # Remove any rows for which there are no input reads:
#   sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
#   kds <- params[1:(nrow(params)-2),]$Mean
#   names(kds) <- rownames(sitesXcounts)
#   num.sf <- nrow(sfXc)

#   kds.s <- params$Mean[1:(nrow(params) - 2)]
#   names(kds.s) <- rownames(params)[1:(nrow(params) - 2)]
#   kd.site <- kds.s[names(kds.s)==site]
#   # Omit the site kd for which the flanking sites are being fit.
#   kds.s <- kds.s[names(kds.s) != site]
#   sitesXcounts <- rbind(sitesXcounts, sfXc)

#   sitesXcounts <- sitesXcounts[rownames(sitesXcounts) != site, ]

#   sitesXcounts <<- sitesXcounts
#   bgs <- rep(params["bg",]$Mean, 5)
#   k.c.stockago <- params["AGO",]$Mean

#   c.I.tots <- Norm(sitesXcounts[,2])*100
#   names(c.I.tots) <- rownames(sitesXcounts)
#   data <- sitesXcounts[,3:7]
#   names(c.I.tots) <- rownames(data)
#   data <<- data
#   x_model <- seq(min(as.numeric(colnames(data))) * 0.7, max(as.numeric(colnames(data))) / 0.7,
#                  length=100)

#   c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
#   colnames(c.totals) <- colnames(x_model)
#   rownames(c.totals) <- rownames(data)

#   c.agos <- sapply(x_model, function(x) {
#       as.numeric(x) * k.c.stockago / 100
#     }
#   )
#   kds.s <<- kds.s
#   flank.kds <<- flank.kds
#   kds <- c(kds.s,flank.kds$Mean)
#   kds <<- kds
#   c.I.tots <<- c.I.tots
#   c.bounds <- as.matrix(
#     sapply(c.agos, function(x) {
#         return(GetBoundRNA(kds, c.I.tots, x))
#       }
#     )
#   )
#   print((c.agos - colSums(c.bounds) )/ c.agos)
#   c.frees <- c.totals - c.bounds
#   c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))

#   if (bgoff == TRUE) {
#     c.all <- c.bounds
#   } else {
#     c.all <- c.bounds + c.bgs
#   }

#   c.final <- data.frame(t(t(c.all) / colSums(c.all)))
#   rownames(c.final) <- rownames(data)
#   x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000/2
#   print(x)
#   y <- c(1,1,1,1,1)
#   print(x)
#   sites.norm <- Norm(c.I.tots)
#   data.norm <- t(t(data)/colSums(data))

#   data.R <- data.norm/(sites.norm)
#   model.R <- c.final/(sites.norm)

#   xmin <- min(x)*0.3
#   xmax <- max(x)*3
#   ymin <- 0.2
#   ymax <- 300
#   yextension <- (ymax/ymin)
#   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   xs <- xs[xs >= xmin & xs <= xmax]
#   xmin <- min(xs)
#   xmax <- max(xs)
#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))
#   ys <- ys[ys >= ymin & ys <= ymax]
#   ymin <- min(ys)
#   ymax <- max(ys)
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))

#   plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.4), ylim=c(ymin, ymax), type="l",
#      col="white", axes=FALSE, ann=FALSE)        
#   # Generate tickmarks for axis.

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0, )
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=xmin, las=2, lwd=0, )
#   axis(2, at=ys, labels=FALSE,
#        pos=xmin)

#   title(main = mirna, font.main=1, line = -2, adj=0.1)
#   title(xlab = "[AGO2-miRNA] (pM)", adj=0.3)
#   title(ylab = "Enrichment")

#   plotlist <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
#                                           "computation/AgoRBNS/",
#                                           "AssignSiteTypes/sites.", mirna,
#                                             "_", plotlist, ".txt"),
#                               stringsAsFactors=FALSE)[,1], "None")

#   legend.names <- rownames(data)[order(kds)]
#   legend.names <- legend.names[which(legend.names %in% plotlist)]
#   ordered_list <- legend.names


#   site_colors <- rep("gray", length(plotlist)-1)
#   flank_colors <- sapply(rownames(flank.kds),GetColorFunction)
#   print(flank_colors)

#   colors_all <- c(site_colors, flank_colors)
#   names(colors_all) <- rownames(data)
#   for (name in rownames(data)) {
#     type = "p"
#     points(x, data.R[name, ], col=colors_all[name], type=type, pch=19,
#            cex=1.2, lwd=3)
#     lines(x_model*k.c.stockago/100*1000/2, model.R[name, ], col=colors_all[name], lwd=2)      
#   }
# }


# PlotBaekKds <- function(mirna, experiment, n_constant, xpos = 20, ypos =20) {
#   # Make plot window.
#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)
#   # Get the kds.
#   kds <- GetSiteKds(mirna, experiment, n_constant, "baek")
#   # Get the list of sites (Redundant with kd names?)
#   site_list <- c(
#     read.table(
#       file = paste0("/lab/bartel1_ata/mcgeary/",
#                     "computation/AgoRBNS/",
#                     "AssignSiteTypes/sites.", mirna,
#                     "_", "baek", ".txt"),
#       stringsAsFactors=FALSE)[,1], "None")
#   # Remove Ago and bg parameters from kd list.
#   kds <- kds[-which(rownames(kds) %in% c("bg", "AGO")),]
#   # Make a data frame where each kd is associated with a sitetype
#   # category.
#   baek_categories <- rownames(kds)
#   baek_categories[which(baek_categories=="5mer-m2.6")] <- "CDNST 1"
#   baek_categories[grep("7mer-A1mm",baek_categories)] <- "CDNST 2"
#   baek_categories[grep("8mer-m2.9", baek_categories)] <- "CDNST 3"
#   baek_categories[grep("8mer-mm", baek_categories)] <- "CDNST 4"
#   kds <- data.frame(Mean = kds$Mean,
#                     Sitetype = baek_categories,
#                     stringsAsFactors = FALSE)
#   # Make dataframe with means for each category.
#   sitetype.means.df = aggregate(
#     kds$Mean, list(kds$Sitetype),function(x) {10^mean(log10(x))})

#   # Get the ordering vector (for colors).
#   site.order = order(sitetype.means.df[,2])
#   # Converts the order to the the per-site positions.

#   # I.E.: If the vector is c(7, 6, 5, 1, etc., then the
#   # new vector will have a 1 at position 7, a 2, at position 6,
#   # a 3 at position 5, a 4 at position 1, etc.)
#   site.rank = sapply(
#     seq(
#       length(site.order)), function(ind) {
#         return(which(site.order == ind))
#     }
#   )
#   par(kPlotParameters)
#   # Initial plot.
#   plot(1, type ="n",
#           axes    = FALSE,
#           ann = FALSE,
#           log = 'x',
#           ylim       = c(0, 22),
#           xlim       = rev(c(0.00003, 3)))
#   # Makes the mirna title.
#   title(main = mirna,
#         line = -2,
#         adj  = 0.1)
#   # Makes the Kd axis title.
#   title(xlab = expression(italic(K)[D]))
  
#   # AXES:
#   ymin=0.0001
#   ymax=3

#   # Makes the axes ticks.
#   # ys refers to the short ticks, of which there are more
#   # yl refers to the long ticks, of which there are fewer.
#   ys <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
#   ys <- ys[which(ys>=ymin & ys <= ymax)]
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#   axis(1, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=0, lwd=0)
#   axis(1, at=ys, labels=FALSE,
#        pos=0)

#   # Plot the beeswarm.
#   beeswarm(Mean ~ Sitetype,
#           data=kds,
#           pch = 21,
#           col = "black",
#           bg = kBaekColors[levels(factor(kds$Sitetype))],
#           method = "swarm",
#           pt.lwd= 1,
#           cex = 1.2,
#           add = TRUE,
#           corral = "random",
#           corralWidth = 0.5,
#           # ylim = c(2,0.00001),
#           at = length(site.rank) -site.rank +1,
#           horizontal = TRUE)

#   # Plots the names of the site types to the right of each category.
#   labelpositions.x <- aggregate(kds$Mean, list(kds$Sitetype), min)[,2]/1.6
#   text(x      = labelpositions.x,
#        y      = length(site.rank)-site.rank +1, 
#        labels = sitetype.means.df[,1],
#        adj    = 0,
#        cex = 0.9,
#        col    = "black")
# }
# SortKdsFile <- function(mirna, sitelist) {
#     # This prints out a new site list that is ordered as per the mean
#     # Kd value from the original sitelist. This allows the Kds to be fit a
#     # second time to make sure that the values are robust to the ordering.
#     kds <- GetSiteKds(mirna, "equilibrium", 5, sitelist)
#     kds <- kds[-which(rownames(kds) %in% c("bg", "AGO", "None")),]

#     write.table(file=paste0("AssignSiteTypes/sites.", mirna, "_", sitelist,
#                             ",ordered.txt"),
#                 x=rownames(kds)[order(kds$Mean)], col.names = FALSE, row.names= FALSE,quote=FALSE)
#     print(kds)
#     print(kds[order(kds$Mean),])
# }


# PlotPositionalKds <- function(experiment, n_constant, sitelist, xpos = 20, ypos = 20) {
#   # Get kds for all site-types of the mirna.

#   positional_sites <- c("8mer", "6mer", "6mer-m8",
#                         "11mer-m3.13",
#                         "11mer-m4.14",
#                         "11mer-m5.15",
#                         "11mer-m6.16",
#                         "11mer-m7.17",
#                         "11mer-m8.18",
#                         "11mer-m9.19",
#                         "11mer-m10.20",
#                         "11mer-m11.21",
#                         "11mer-m12.22",
#                         "11mer-m13.23",
#                         "None")
#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)
#   kds <- GetSiteKds("miR-1", experiment, n_constant, sitelist)    
#   kds <- kds[positional_sites,]
#   print(kds)
#   par(kPlotParameters)

#   # Make plot with miR-1 data:
#   ind_p <- c(1,2,3,nrow(kds))
#   y_p <- (nrow(kds) - seq(nrow(kds)) + 1)[ind_p]
#   ind_l <- seq(nrow(kds))[-ind_p]
#   y_l <- (nrow(kds) - seq(nrow(kds)) + 1)[ind_l]
#   plot(1, type = "n",
#         col = "white",
#         axes    = FALSE,
#         log = 'x',
#         ylim       = c(0, 22),
#         xlim       = rev(c(0.00003, 5)))

#   title(xlab = expression(italic(K)[D]))
#   ymin=0.0001
#   ymax=5
#   ys <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
#   ys <- ys[which(ys>=ymin & ys <= ymax)]
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#   axis(1, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=0)
#   axis(1, at=ys, labels=FALSE,
#        pos=0)

#   for (mirna in c("miR-1","let-7a", "miR-124", "lsy-6", "miR-155")) {
#     kds <- GetSiteKds(mirna, experiment, n_constant, sitelist)    
#   kds <- kds[positional_sites,]
#   points(kds$Mean[ind_p], y_p,
#         col = kMirnaColors[mirna],
#         pch=19)
#   lines(kds$Mean[ind_l], y_l,
#         col = kMirnaColors[mirna],
#         lwd = 2,
#         type = "o")
#   lines(kds$Lower_CI[ind_l], y_l,
#         col = kMirnaColors[mirna],
#         lwd = 1,
#         lty=2,
#         cex = 1)
#   lines(kds$Upper_CI[ind_l], y_l,
#         col = kMirnaColors[mirna],
#         lwd = 1,
#         lty=2,
#         cex = 1)
#   }
#   text(6.5,
#        nrow(kds) - seq(nrow(kds)) + 1,
#        labels=rownames(kds),
#        cex=0.9,
#        adj=0,
#        col= "black")
#   legend(x = 10^-3,
#          y = 8,
#          legend = mirnas.all,
#          col = kMirnaColors,
#          bty="n",
#          pch = 19)
# }

# PlotPositionalKdsMiR7 <- function(n_constant, sitelist, xpos = 20, ypos = 20) {
#   # Get kds for all site-types of the mirna.

#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)
#   kds <- GetSiteKds("miR-7-22nt", "equilibrium_nb", n_constant, sitelist)    
#   print(kds)
#   par(kPlotParameters)

#   # Make plot with miR-1 data:
#   ind_p <- c(1,2,3,nrow(kds))
#   y_p <- (nrow(kds) - seq(nrow(kds)) + 1)[ind_p]
#   ind_l <- seq(nrow(kds))[-ind_p]
#   y_l <- (nrow(kds) - seq(nrow(kds)) + 1)[ind_l]
#   plot(1, type = "n",
#         col = "white",
#         axes    = FALSE,
#         log = 'x',
#         ylim       = c(0, 22),
#         xlim       = rev(c(0.00003, 5)))

#   title(xlab = expression(italic(K)[D]))
#   ymin=0.0001
#   ymax=5
#   ys <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
#   ys <- ys[which(ys>=ymin & ys <= ymax)]
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#   axis(1, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=0)
#   axis(1, at=ys, labels=FALSE,
#        pos=0)
#   positional_sites <- rownames(kds)[1:(nrow(kds) - 2)]
#   kMirnaColors <- topo.colors(8)[2:5]
#   names(kMirnaColors) <- c("miR-7-22nt","miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
#   mirnas.all <- names(kMirnaColors)
#   for (mirna in c("miR-7-22nt","miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
#     kds <- GetSiteKds(mirna, "equilibrium_nb", n_constant, sitelist)    
#   kds <- kds[positional_sites,]
#   print(kds)
#   points(kds$Mean[ind_p], y_p,
#         col = kMirnaColors[mirna],
#         pch=19)
#   lines(kds$Mean[ind_l], y_l,
#         col = kMirnaColors[mirna],
#         lwd = 2,
#         type = "o")
#   lines(kds$Lower_CI[ind_l], y_l,
#         col = kMirnaColors[mirna],
#         lwd = 1,
#         lty=2,
#         cex = 1)
#   lines(kds$Upper_CI[ind_l], y_l,
#         col = kMirnaColors[mirna],
#         lwd = 1,
#         lty=2,
#         cex = 1)
#   }
#   text(6.5,
#        nrow(kds) - seq(nrow(kds)) + 1,
#        labels=rownames(kds),
#        cex=0.9,
#        adj=0,
#        col= "black")
#   legend(x = 10^-3,
#          y = 8,
#          legend = mirnas.all,
#          col = kMirnaColors,
#          bty="n",
#          pch = 19)
# }


# PlotSiteFlankKds <- function(mirna, experiment, n_constant, sitelist,
#                                 plotlist, xpos = 20, ypos = 20) {

#   if (length(plotlist) == 0) {
#     site_list <- rownames(data)
#   } else if (class(plotlist) == "character") {
#     site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
#                                             "computation/AgoRBNS/",
#                                             "AssignSiteTypes/sites.", mirna,
#                                             "_", plotlist, ".txt"),
#                               stringsAsFactors=FALSE)[,1], "None")
#   } else {
#     site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
#   }
#   num.sites <- length(site_list)
#   # Get kds for all site-types of the mirna.
#   params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
#   site.kds <- params$Mean[1 : (nrow(params) - 2)]
#   names(site.kds) <- rownames(params)[1 : (nrow(params) - 2)]
#   # # Use the 8mer flanks to just get the flank strings.
#   flanks.temp <- GetFlankKds(mirna, experiment, n_constant, sitelist, "8mer")
#   # Get the flank string names
#   flanks <- rownames(flanks.temp)

#   # Pre-allocate the matrix with the flanking kds.
#   flank.kds <- matrix(NaN,nrow=length(flanks),ncol=length(site_list)-1)
#   # Name the rows and columns.
#   rownames(flank.kds) <- flanks
#   colnames(flank.kds) <- site_list[-length(site_list)]
#   # Remove the teporary flank data structure.
#   rm(flanks.temp)

#   # Fill in the flanking kd matrix for each sitetype and flanking dinucleotide
#   # combination that exists.
#   for (site in site_list[1:(length(site_list)-1)]) {
#     flank.kds.site <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
#     for (flank in rownames(flank.kds.site)) {
#       flank.kds[flank,site] <- flank.kds.site[flank,]$Mean
#     }
#   }
#   # Removes site types for which there are literally zero flank Kds.
#   site.kds <- site.kds[which(colSums(is.na(flank.kds))!=256)]
#   site.kds <- site.kds[site_list[-length(site_list)]]
#   flank.kds.trim <- flank.kds[,site_list[-length(site_list)]]
#   flank.kds.trim <- flank.kds.trim[,order(site.kds)]
#   kds.flanks <- c(flank.kds.trim)
#   data.sites <- rep(colnames(flank.kds.trim), each=nrow(flank.kds.trim))
#   data.ranks <- rep(num.sites-seq(num.sites-1), each=nrow(flank.kds.trim))
#   data.colors <- rep(sapply(rownames(flank.kds), GetColorFunction, alpha=1),
#                      ncol(flank.kds.trim))
#   flanks.df <- data.frame(kds=log10(kds.flanks), rank=as.numeric(data.ranks), sites=data.sites,
#                           cols=data.colors,
#                            stringsAsFactors=FALSE)

#   print(unique(flanks.df$site))
#   flanks.df <<- flanks.df
#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)

#   par(kPlotParameters)
#   boxplot(kds ~ rank,
#           data       = flanks.df,
#           axes       = FALSE,
#           horizontal = TRUE,
#           outline    = FALSE,
#           xlim       = c(0, 22),
#           ylim       = log10(rev(c(0.00003, 10))))
#   print(c(0,num.sites+5))
#   title(main = mirna,
#         line = -2,
#         adj  = 0.1)
#   title(xlab = expression(italic(K)[D]))
#     beeswarm(kds ~ rank,
#              data       = flanks.df,
#              add        = TRUE,
#              method     = "swarm",
#              corral     = "random",
#              corralWidth = 0.5,
#              pch        = 1,
#              cex = 0.8,
#              horizontal = TRUE,  
#              pwcol = cols)




#   ymin=0.0001
#   ymax=10
#   ys <- c(sapply(seq(floor(min(flanks.df$kds, na.rm=TRUE)),
#                      ceiling(max(flanks.df$kds, na.rm=TRUE))), function(x) seq(10)*10^x))
#   ys <- ys[ys >= ymin & ys <= ymax]
#   ymin <- min(ys)
#   ymax <- max(ys)
#   yl <- seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#   axis(1, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=0, lwd=0)
#   axis(1, at=log10(ys), labels=FALSE,
#        pos=0)
#   # points(rep(-3.4, num.sites-1),num.sites - seq(num.sites-1),col= kSiteColors[names(site.kds),][order(site.kds)],cex=1.5,pch=19)

#   flank_mins <- apply(flank.kds.trim,2,function(col) {min(col[!is.na(col)])})
#     text(log10(flank_mins) - 0.2,
#          num.sites - seq(num.sites-1),
#          labels=ConvertTtoU(unique(flanks.df$sites)),
#          adj=0, col= "black", cex = 0.8)

#   arrows(log10(max(flank.kds[,"8mer"])), 16.8,
#          log10(min(flank.kds[,"8mer"])), 16.8,
#          length=0.05, lty=1, angle=90, code=3)
#   text(-2.1,
#        17.3, "Range of 8mer binding affinity across dinucleotide context", cex = 0.6)

#   arrows(log10(site.kds["6mer"]), 18,
#         log10(site.kds["8mer"]), 18,
#         length=0.05, lty=1, angle=90, code=3)
#   text(-2.1,
#        18.5, "Difference in average binding betwen 6mer and 8mer site", cex = 0.6)

# }

# PlotSiteFlankKdsScatter <- function(mirna, experiment, n_constant_1, n_constant_2, sitelist,
#                                 site, xpos = 20, ypos = 20) {


#   flanks.kds.1 <- GetFlankKds(mirna, experiment, n_constant_1, sitelist, "8mer")
#   flanks.kds.2 <- GetFlankKds(mirna, experiment, n_constant_2, sitelist, "8mer")

#   # Get the flank string names


#   dev.new(xpos = xpos, ypos = ypos, height = 5, width = 5)
#     par(kPlotParameters)

#   plot(flanks.kds.1$Mean,
#     flanks.kds.2$Mean,
#     xlim = c(0.0001, 5),
#     ylim = c(0.0001, 5),
#       axes = FALSE,
#       col=sapply(rownames(flanks.kds.1), GetColorFunction, alpha=1),
#       ann= FALSE,
#       log ='xy')
#   title(main = mirna,
#         line = -2,
#         adj  = 0.1)
#   title(xlab = expression(italic(K)[D]))


#   ymin=0.0001
#   ymax=5


#   ys <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
#   ys <- ys[which(ys>=ymin & ys <= ymax)]
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#   print(yl)
#   print(ys)
#   axis(1, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin)
#   axis(1, at=ys, labels=FALSE,
#        pos=ymin)

#   # points(rep(-3.4, num.sites-1),num.sites - seq(num.sites-1),col= kSiteColors[names(site.kds),][order(site.kds)],cex=1.5,pch=19)


# }

# GetCanonicalSiteFlankingKdMatrix <- function(experiment, n_constant, sitelist) {

#   flanks.temp <- GetFlankKds("miR-1", experiment, n_constant, sitelist, "8mer")
#   flank.strings <- rownames(flanks.temp)

#   data.matrix <- matrix(NA, ncol = 258, nrow =30)
#   colnames(data.matrix) <- c("miRNA", "site", flank.strings)
#   row <- 1
#   for (mirna in mirnas.all) {
#     for (site in canonical.sites) {
#       flanks <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
#       data.matrix[row, 1] <- mirna
#       data.matrix[row, 2] <- site
#       row_mean <- mean(log10(as.numeric(flanks$Mean)))
#       for (flank in rownames(flanks)) {
#         data.matrix[row,which(colnames(data.matrix)==flank)] <- log10(as.numeric(flanks[flank,]$Mean)) - row_mean
#       }
#       row <- row + 1
#     }
#   }
#   return(data.matrix)
# }

# PlotCanonicalSiteFlankingKdHeatmap <- function(experiment, n_constant, sitelist, xpos = 20, ypos = 20) {
#   data.matrix <- GetCanonicalSiteFlankingKdMatrix(experiment, n_constant, sitelist)
#   data.flanks <- data.frame(data.matrix[,3:258],stringsAsFactors=FALSE)
#   rownames(data.flanks) <- apply(data.matrix,1,function(row) {
#     paste0(c(row[1],row[2]), collapse=" | ")
#     })
#   dev.new(xpos = xpos, ypos = ypos, width = 14, height = 5)
#   heatmap.2(data.matrix(data.flanks),
#     col=rainbow(100,s=0.8,v=1,start=0.40,end=.1),
#     Rowv = "NULL",
#     dendrogram="column",
#     trace = "none",
#     scale = "none",
#     cexRow = 0.9,
#     margins = c(3.5, 8),
#     na.color="gray",
#     lhei=c(3,10),
#     lwid = c(3, 9),
#     key=TRUE)

# }

# GenerateLMForFlankingKds <- function(experiment, n_constant, sitelist) {
#   data.matrix <- GetCanonicalSiteFlankingKdMatrix(experiment, n_constant, sitelist)
#   print(data.matrix[1:5, 1:5])
#   print(c(data.matrix[3:5,3:5]))
# }

# GetPairingFlankData <- function(mirna, experiment, condition, n_constant,
#                                 sitelist, site, mir.start, mir.stop,
#                                 absolute=TRUE, noconstant=FALSE) {
#   if (absolute == TRUE) {
#     folder = "/structural_analysis_PAPER_realfinal/"
#   } else {
#     folder = "/structural_analysis_PAPER_final/"
#   }
#   file = paste0(condition, "_", n_constant, "_", mir.start, "-", mir.stop,
#                 collapse="")
#   if (noconstant == TRUE) {
#     file = paste0(file, "_noconstant", collapse="")
#   }
#   data <- read.table(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
#                             experiment, folder, site, "/", file, ".txt"),
#                      header=TRUE)
#   data <- data.frame(data)
#   return(data)
# }


# MakeFlankScatterPlotSimple <- function(mirna, experiment, condition, n_constant, sitelist, site, mir.start, mir.stop) {
#   data <- data.frame(GetPairingFlankData(mirna, experiment, condition, n_constant, sitelist, site, mir.start, mir.stop))

#   print(data[1:10,])
#   print(data$Flank)
#   attach(data)
#   par(kPlotParameters)
#   plot(data$plfold^(1/(mir.stop - mir.start + 1)),
#        data$accessibility,
#        col=sapply(unlist(lapply(data$Flank,as.character)), GetColorFunction, alpha=0.1))
#   detach(data)
# }

# MakeFlankScatterPlot <- function(mirna, experiment, condition, n_constant, sitelist, site, mir.start, mir.stop) {
#   data <- data.frame(GetPairingFlankData(mirna, experiment, condition, n_constant, sitelist, site, mir.start, mir.stop))

#   print(data[1:10,])
#   print(data$Flank)
#   attach(data)
#   par(kPlotParameters)
#   plot(by(data,Flank,function(x){mean(x$plfold^(1/(mir.stop - mir.start + 1)))}),
#        by(data,Flank,function(x){mean(x$accessibility)}),
#        col=sapply(unique(unlist(lapply(data$Flank,as.character))), GetColorFunction, alpha=1),
#        log='xy',
#        xlim = c(0.03,1),
#        ylim = c(0.03,1))
#   detach(data)
# }

# PlotStructECDF <- function(mirna, experiment, n_constant, sitelist, site,
#                            mir.start, mir.stop, condition = FALSE, score="plfold", flank=FALSE,
#                            absolute=TRUE, logparam="",noconstant=FALSE,
#                            xpos=20, ypos = 20, xlims=c(0.2, 1),
#                            ylims=c(0, 1)) {
#   dev.new(xpos=xpos, ypos=ypos, height=5, width=5)

#   xs     <- seq(xlims[1], xlims[2], by=0.2)
#   x_ecdf <- seq(xlims[1], xlims[2], by=0.002)

#   par(kPlotParameters)
#   par(lwd=1.5)
#   plot(1, type = "n",
#        xlim = xlims,
#        ylim = ylims,
#        log  = logparam,
#        axes = FALSE,
#        ann  = FALSE)

#   axis(1, at=xs, pos=ylims[1])
#   axis(2, at=seq(ylims[1], ylims[2], by=0.2), pos=xlims[1])
#   if (condition != FALSE){
#     conditions = c("I_combined", condition)
#   } else {
#     conditions = c("I_combined", "0", "0.4", "1.26", "4", "12.6", "40")
#   }
#   sapply(conditions,
#          function(condition) {
#     data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
#                                            n_constant, sitelist, site,
#                                            mir.start, mir.stop,
#                                            absolute=absolute,
#                                            noconstant=noconstant))
#     if (flank != FALSE) {
#       data <- subset(data, Flank==flank)
#     }
#     data.intermediate <- data[[score]]
#     ecdf.d <- ecdf(data.intermediate)
#     lines(x_ecdf, ecdf.d(x_ecdf), col=kEquilCondColors[condition])
#   })

#   title(xlab = "Per-nucleotide accessibility across from miRNA nt 1-15")
#   title(ylab = "Cumulative Fraction", line = 1)

#   legend.names <- c("Input library",
#                     "0% AG02-miR-1",
#                     "0.4% AG02-miR-1",
#                     "1.26% AG02-miR-1",
#                     "4% AGO2-miR-1",
#                     "12.6% AGO2-miR-1",
#                     "40% AGO2-miR-1")
#   names(legend.names) <- c("I_combined", "0", "0.4", "1.26", "4", "12.6", "40")

#   legend(x=0.225, y=1, cex=0.9, legend=legend.names[conditions], lwd=1, col=kEquilCondColors[conditions],
#          bty="n", ncol=1)
# }



# PlotStructureVsFlankingKds <- function(mirna, experiment, condition, n_constant,
#                                       sitelist, site, mir.start, mir.stop,
#                                       xpos = 20, ypos = 20, absolute=TRUE,
#                                       noconstant=TRUE) {
#   dev.new(xpos = xpos, ypos = ypos, width = 5, height = 5)
#   par(kPlotParameters)
#   # Window size for normalizing the pl_fold with the accessibility:
#   win = 1/(mir.stop - mir.start + 1)

#   # Get the flanking kds:
#   flank.kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
#   kds <- flank.kds$Mean
#   names(kds) <- sapply(rownames(flank.kds), function(name) {
#     paste0(unlist(strsplit(name, split = ""))[-3], collapse = "")
#     })

#   plot_subfunction <- function(condition,factor, ylimit, logparam=FALSE) {
#         data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
#                            n_constant, sitelist, site, mir.start, mir.stop,
#                            absolute=absolute, noconstant=noconstant))
#     print(data[1 : 10, ])

#     p_access <- c(by(data, data$Flank, function(x) {
#         mean(x[[factor]])
#     }))
#     xmin <- 0.0003
#     xmax <- 0.2

#     if (condition == "I_combined") {
#       fit <<- lm(log10(p_access)~log10(kds))

#     }
#     plot(1, type = "n",
#          log = 'xy',
#          xlim = rev(c(xmin , xmax)),
#          ylim = ylimit,
#          axes = FALSE,
#          ann = FALSE)
#     ymin=ylimit[1]
#     ymax=ylimit[2]
#     m <- fit$coefficients[2]
#     b <- fit$coefficients[1]
#     x_line <- 10^seq(log10(xmin), log10(xmax),length = 20)
#     y_line <- 10^(m*log10(x_line) + b)

#     lines(x_line, y_line, lty = 2,lwd = 0.5)

#     points(kds,
#          p_access[names(kds)],
#          col=sapply(names(kds), GetColorFunction, alpha=0.1),
#          cex = 0.8)
#     xs <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
#   xs <- xs[which(xs>=xmin & xs <= xmax)]
#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))

#   ys <- seq(ylimit[1], ylimit[2], by = 0.1)
#   yl <- seq(ylimit[1], ylimit[2], by = 0.2)
#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin, lwd=1)

#   axis(2, at=yl,
#        pos=xmax, lwd=0)
#   axis(2, at=ys, labels=FALSE,
#        pos=xmax, lwd=1)
#   title(xlab = expression(italic(K)[D]))
#   title(ylab = "Average accessibility across from miRNA nt 1-15", line=1.5)


#     if (condition == "I_combined") {
#       condition <- "I"
#     }
#     text(x=0.1, y=0.95, condition)
#       cor_text <- format(round(cor(-log(kds),log(p_access[names(kds)]),
#                                 use = "pairwise.complete.obs")^2,2), nsmall = 2)     

#     text(x = 0.001, y = 0.95, bquote(italic(r)^2 ==  .(cor_text)))

#   }

#   plot_subfunction("I_combined", factor = "plfold", ylimit = c(0.4,1), logparam=FALSE)

# }

# MakeAccessibilityvsInputPlot <- function(mirna, experiment, condition, n_constant,
#                                       sitelist, site, mir.start, mir.stop,pl=TRUE,
#                                       noconstant=FALSE, absolute=TRUE) {
#   dev.new(xpos = 20, ypos = 220, height = 5, width = 6)
#   par(kPlotParameters)
#   # Window size for normalizing the pl_fold with the accessibility:
#   win = 1/(mir.stop - mir.start + 1)

#   # Get the flanking kds:
#   flank.kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
#   kds <- flank.kds$Mean
#   names(kds) <- sapply(rownames(flank.kds), function(name) {
#     paste0(unlist(strsplit(name, split = ""))[-3], collapse = "")
#     })
#     subfunction <- function(condition_temp) {
#     data <- data.frame(GetPairingFlankData(mirna, experiment, condition_temp,
#                        n_constant, sitelist, site, mir.start, mir.stop,
#                        noconstant=noconstant, absolute=absolute))
#     attach(data)

#     p_access<- c(by(data, Flank, function(x) {
#           10^(mean(log(x$plfold^win)))
#       }
#     ))
#     detach(data)

#     return(p_access)
#   }
#     p_I <- subfunction("I")
#     p_condition <- subfunction(condition)
#     plot(p_I,
#          p_condition[names(p_I)],
#          col=sapply(names(p_I), GetColorFunction, alpha=0.1),
#          log = 'xy',
#          xlim = c(0.1, 1),
#          ylim = c(0.1, 1))
  

 

# }

# PlotPlByAUScoreBin <- function(mirna, experiment, condition, n_constant,
#                                    sitelist, site, mir.start, mir.stop,
#                                    xpos=20, ylimit = c(0, 1), xlimit = c(0, 1),
#                                    num.bins=10, inside.range = c(0.3, 0.8),
#                                    ypos=20, AU_score="AU_win",
#                                    p_score="plfold", noconstant=FALSE,
#                                    absolute=TRUE) 
# {
#   dev.new(xpos=xpos, ypos=ypos, width=5, height=5)
#   par(kPlotParameters)

#   xs <- seq(xlimit[1], xlimit[2], by = 0.2)
#   ys <- seq(ylimit[1], ylimit[2], by = 0.1)
#   plot(1, type = "n",
#        xlim = xlimit,
#        ylim = ylimit,
#        ann = FALSE,
#        axes = FALSE)

#   bin.limits <- c(0, seq(inside.range[1], inside.range[2], length = num.bins), 1)
#   print(bin.limits)
#   b.l <- length(bin.limits)
#   x <- (bin.limits[2 : b.l] + bin.limits[1 : (b.l - 1)]) / 2
#   print(x)
#   # Window size for normalizing the pl_fold with the accessibility:
#   subfunction <- function(condition, factor, cond.color = "black") {
#     data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
#                        n_constant, sitelist, site, mir.start, mir.stop,
#                        noconstant=noconstant, absolute=absolute))

#     data$plfold_bins = cut(data[[AU_score]], bin.limits,labels=x)

#     agg <- aggregate(. ~ plfold_bins, data, function(x) {
#       c(mean = mean(x),
#         se = sd(x)/sqrt(length(x)-1))
#       })

#     x <- as.numeric(as.character(agg$plfold_bins))
#     y <- agg[[p_score]][,1]
#     y.e <- agg[[p_score]][,2]
#     segments(x0=xlimit[1], y0=mean(data[[p_score]]), x1=xlimit[2], y1=mean(data[[p_score]]),col=cond.color,lty=2)

#     points(x, y,
#          type = "o",
#          col=cond.color)
#     arrows(x, y - y.e, x, y + y.e,
#            col=cond.color, length=0.15, angle=90, code=3)

#   }




#     subfunction("I_combined", AU_score)
#     subfunction(condition, AU_score,cond.color="forestgreen")
#     axis(1, xs, pos = ylimit[1])
#     axis(2, ys, pos = xlimit[1])
#     title(xlab = "Mean AU context-score across from miRNA nt 1-15")
#     title(ylab = "Per-nucleotide accessibility across from miRNA nt 1-15", line=1)
#     legend(x=0.025, y = 1, bty = "n", legend = c("Input library", "4% AGO2-miR-1"), col = c("black","forestgreen"), lwd=1)

# }

# PlotAUScoreByPlBin <- function(mirna, experiment, condition, n_constant,
#                                  sitelist, site, mir.start, mir.stop, 
#                                  noconstant=FALSE, absolute=TRUE, xpos=20,
#                                  ypos=20, xlimit=c(0,1), ylimit=c(0,1),num.bins=10,
#                                  inside.range = c(0.3, 0.8),
#                                  AU_score="AU_win", p_score="plfold")
# {
#   dev.new(xpos=xpos, ypos=ypos, width=5, height=5)
#   par(kPlotParameters)

#   xs <- seq(xlimit[1], xlimit[2], by = 0.2)

#   ys <- seq(ylimit[1], ylimit[2], by = 0.05)
#   plot(1, type = "n",
#        xlim = xlimit,
#        ylim = ylimit,
#        ann = FALSE,
#        axes = FALSE)

#   bin.limits <- c(0, seq(inside.range[1], inside.range[2], length = num.bins), 1)
#   print(bin.limits)
#   b.l <- length(bin.limits)
#   x <- (bin.limits[2 : b.l] + bin.limits[1 : (b.l - 1)]) / 2
#   print(x)
#   # Window size for normalizing the pl_fold with the accessibility:
#   subfunction <- function(condition, factor, cond.color = "black") {
#     data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
#                        n_constant, sitelist, site, mir.start, mir.stop,
#                        noconstant=noconstant, absolute=absolute))

#     data$plfold_bins = cut(data[[p_score]], bin.limits,labels=x)

#     agg <- aggregate(. ~ plfold_bins, data, function(x) {
#       c(mean = mean(x),
#         se = sd(x)/sqrt(length(x)-1))
#       })

#     x <- as.numeric(as.character(agg$plfold_bins))
#     y <- agg[[AU_score]][,1]
#     y.e <- agg[[AU_score]][,2]
#     segments(x0=xlimit[1], y0=mean(data[[AU_score]]), x1=xlimit[2], y1=mean(data[[AU_score]]),col=cond.color,lty=2)

#     points(x, y,
#          type = "o",
#          col=cond.color)
#     arrows(x, y - y.e, x, y + y.e,
#            col=cond.color, length=0.15, angle=90, code=3)


#   }




#     subfunction("I_combined", AU_score)
#     subfunction(condition, AU_score,cond.color="forestgreen")
#     axis(1, xs, pos = ylimit[1])
#     axis(2, ys, pos = xlimit[1])
#     title(xlab = "Per-nucleotide accessibility across from miRNA nt 1-15")
#     title(ylab = "Mean AU context-score across from miRNA nt 1-15", line=1)
#     legend(x=0.025, y = 0.6, bty = "n", legend = c("Input library", "4% AGO2-miR-1"), col = c("black","forestgreen"), lwd=1)

# }

# MakeDoubleSubsetMatrix <- function(mirna, experiment, condition, n_constant,
#                                  sitelist, site, mir.start, mir.stop, 
#                                  noconstant=FALSE, absolute=TRUE, xpos=20,
#                                  ypos=20, xlimit=c(0,1), ylimit=c(0,1),num.bins.AU=8,
#                                  num.bins.pl=6,
#                                  inside.range.AU = c(0.3, 0.7),
#                                  inside.range.pl = c(0.3,0.7),
#                                  AU_score="AU_cs", p_score="plfold")
# {
#   dev.new(xpos=xpos, ypos=ypos, width=10, height=10)

#   par(kPlotParameters)
#   par(mfrow=c(2,2))
#   xs <- seq(xlimit[1], xlimit[2], by = 0.2)

#   ys <- seq(ylimit[1], ylimit[2], by = 0.05)
#   # plot(1, type = "n",
#   #      xlim = xlimit,
#   #      ylim = ylimit,
#   #      ann = FALSE,
#   #      axes = FALSE)

#   bin.limits.AU <- c(0, seq(inside.range.AU[1], inside.range.AU[2], length = num.bins.AU), 1)
#   bin.limits.pl <- c(0, seq(inside.range.pl[1], inside.range.pl[2], length = num.bins.pl), 1)

#   print(bin.limits)
#   b.l.pl <- length(bin.limits.pl)
#   x.pl <- (bin.limits.pl[2 : b.l.pl] + bin.limits.pl[1 : (b.l.pl - 1)]) / 2

#   b.l.AU <- length(bin.limits.AU)
#   x.AU <- (bin.limits.AU[2 : b.l.AU] + bin.limits.AU[1 : (b.l.AU - 1)]) / 2


#   # Window size for normalizing the pl_fold with the accessibility:
#   subfunction <- function(condition, factor, cond.color = "black") {
#     data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
#                        n_constant, sitelist, site, mir.start, mir.stop,
#                        noconstant=noconstant, absolute=absolute))

#     data$plfold_bins = cut(data[[p_score]], bin.limits.pl,labels=x.pl)
#     data$AU_bins = cut(data[[AU_score]], bin.limits.AU,labels=x.AU)

#     output = matrix(0,nrow=length(levels(data$plfold_bins)),
#                       ncol=length(levels(data$AU_bins)))
#     rownames(output) <- as.character(levels(data$plfold_bins))
#     colnames(output) <- as.character(levels(data$AU_bins))
#     for (pl_bin in levels(data$plfold_bins)) {
#       for (AU_bin in levels(data$AU_bins)) {
#         count <- nrow(subset(data, (plfold_bins == pl_bin & AU_bins == AU_bin)))+0.1
#         output[as.character(pl_bin), as.character(AU_bin)] <- count
#       }
#     }
#     output <- output/sum(output)
#     # x <- as.numeric(as.character(agg$plfold_bins))
#     # y <- agg[[AU_score]][,1]
#     # y.e <- agg[[AU_score]][,2]
#     # segments(x0=xlimit[1], y0=mean(data[[AU_score]]), x1=xlimit[2], y1=mean(data[[AU_score]]),col=cond.color,lty=2)

#     # points(x, y,
#     #      type = "o",
#     #      col=cond.color)
#     # arrows(x, y - y.e, x, y + y.e,
#     #        col=cond.color, length=0.15, angle=90, code=3)

#     return(output)
#   }




#     matrix_I <- subfunction("I_combined", AU_score)
#     print(matrix_I)
#     matrix_A <- subfunction(condition, AU_score,cond.color="forestgreen")
#     print(matrix_A)
#     final <- matrix_A / matrix_I
#     rownames(final) <- round(as.numeric(rownames(matrix_A)), 2)
#     colnames(final) <- round(as.numeric(colnames(matrix_A)), 2)
#     # heatmap.2(final,Rowv=FALSE, Colv=FALSE,
#     #       trace = "none",
#     # cexRow = 0.9,
#     # margins = c(3.5, 3),
#     # na.color="gray",
#     # lhei=c(1,10),
#     # lwid = c(1, 9),
#     # key=FALSE)
#     plot(colnames(final), final[1,], col = "white",xlim = c(0, 1), ylim=c(0.01, 5),type="l")
#         colramp = rainbow(nrow(final),start=0,end=0.5)

#     for (row in 1:nrow(final)){
#       lines(colnames(final), final[row,],col=colramp[row])
#     }

#     plot(rownames(final), final[,1],col="white", xlim = c(0, 1), ylim=c(0.01, 5),type="l")
#     colramp = rainbow(ncol(final),start=0,end=0.5)
#     for (column in 1:ncol(final)){
#       lines(rownames(final), final[,column],col=colramp[column])
#     }

#     plot(colnames(final), final[1,], col = "white",xlim = c(0, 1),log='y', ylim=c(0.01, 5),type="l")
#         colramp = rainbow(nrow(final),start=0,end=0.5)

#     for (row in 1:nrow(final)){
#       lines(colnames(final), final[row,],col=colramp[row])
#     }

#     plot(rownames(final), final[,1],col="white", xlim = c(0, 1), log='y',ylim=c(0.01, 5),type="l")
#     colramp = rainbow(ncol(final),start=0,end=0.5)
#     for (column in 1:ncol(final)){
#       lines(rownames(final), final[,column],col=colramp[column])
#     }
#     dev.new()
#     heatmap.2(final,Rowv=FALSE, Colv=FALSE,
#           trace = "none",
#     cexRow = 0.9,
#     margins = c(3.5, 3),
#     na.color="gray",
#     lhei=c(1,10),
#     lwid = c(1, 9),
#     key=FALSE)

# }



# PlotPlSampleFlanks <- function(mirna, experiment, condition, n_constant,
#                                  sitelist, site, mir.start, mir.stop, noconstant=FALSE,
#                                  absolute=TRUE, xpos=20, depth=5000,
#                                  ypos=20, p_score="plfold")
# {
#   dev.new(xpos=xpos, ypos=ypos , height=5, width=8)
#   par(kPlotParameters)
#   layout(matrix(c(1, 2), 1, 2, byrow=TRUE), widths = c(1.5, 1))
#   # Window size for normalizing the pl_fold with the accessibility:
#   data.I <- data.frame(GetPairingFlankData(mirna, experiment, "I_combined",
#                        n_constant, sitelist, site, mir.start, mir.stop,
#                        noconstant=noconstant, absolute=absolute))

#   data.A <- data.frame(GetPairingFlankData(mirna, experiment, condition,
#                        n_constant, sitelist, site, mir.start, mir.stop,
#                        noconstant=noconstant, absolute=absolute))
#   NumFlanks <- function(data) {
#     flank_abundances <- aggregate(. ~ Flank, data, function(x) {
#       length(x)
#     })
#     flank_num <- flank_abundances[,2]
#     names(flank_num) <- flank_abundances$Flank
#     return(flank_num)
#   }
#   sample_flank_probabilities <- function(dist1,dist2) {
#     weights <- NumFlanks(dist2)/NumFlanks(dist1)
#     inds <- sample(1:nrow(dist1),
#                      replace = TRUE,
#                      size = nrow(dist2),
#                      prob = weights[dist1$Flank])
#     return(inds)
#   }
#   print(model[site,])

#   sample_distribution <- function(dist1, dist2, score) {

#     # Make vector to sample from dists 1 and 2
#     input <- dist1[[score]]
#     len <- length(input)
#     len2 <- nrow(dist2)
#     # Make the target mean:
#     target <- mean(log10(dist2[[score]]^15))
#     target2 <- sd(log10(dist2[[score]]^15))
#     # target <- mean(dist2[[score]])
#     # Define Cost function:
#     size_sample <- as.integer(model[site,as.character(condition)]*len)
#     tick <- 1
#     costfunction_plfold <- function(pars) {
#       tick <<- tick + 1
#       inds <- c()
#       while (length(inds) < len/10){
#         inds <- c(inds,
#                   sample_int_rej(len,
#                                size = as.integer(Logistic(pars[1],1)*len),
#                                prob = input^pars[2] + 10^pars[3]))

#       }
#       residual <- (
#         (mean(log10(input[inds]^15)) - target)^2
#         )
#       return(residual)
#     }
#     pars <- optim(c(0,1,-10),costfunction_plfold)$par
#     inds <- c()
#     while (length(inds) < len2) {
#       inds <- c(inds,
#                 sample_int_rej(len,
#                              size = as.integer(Logistic(pars[1],1)*len),
#                              prob = input^pars[2] + 10^pars[3]))
#       }
#     return(inds)
#     }
#     inds_AU <- sample_flank_probabilities(data.I, data.A)

#   #Plot the probabilities:
#   IndsEffect <- function(condition) {
#     print(condition)
#     data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
#                        n_constant, sitelist, site, mir.start, mir.stop,
#                        noconstant=noconstant, absolute=absolute))
#     inds <- sample_flank_probabilities(data.I, data)
#     return((mean(log10(data.I[[p_score]][inds]^15)) -
#             mean(log10(data.I[[p_score]]^15))
#             )/(mean(log10(data[[p_score]]^15)) -
#                mean(log10(data.I[[p_score]]^15))))
#     }
#         xs = 10^seq(-8,0,by = .2)
#         xl = 10^seq(-8,0,by = 2)

#   plot(1, type = "n",
#          xlim = c(min(xs), 1),
#          log = 'x',
#          ylim = c(0, 1),
#          axes=FALSE,
#          ann=FALSE)
#   x_ecdf <- 10^seq(-8,0,by = 0.01)
#   lines(x_ecdf,ecdf(data.I[[p_score]]^15)(x_ecdf),lwd = 1)
#   lines(x_ecdf,ecdf(data.A[[p_score]]^15)(x_ecdf), col = kEquilCondColors[as.character(condition)],lwd = 2)
#   lines(x_ecdf,ecdf(data.I[[p_score]][inds_AU]^15)(x_ecdf), col = kEquilCondColors[as.character(condition)],lty=2)
  
#   print(c(mean(data.I[[p_score]]^15),
#            mean(data.I[[p_score]][inds_AU]^15),
#            mean(data.A[[p_score]]^15)))
#   print(c(ecdf(data.I[[p_score]]^15)(mean(data.I[[p_score]]^15)),
#            mean(data.I[[p_score]][inds_AU]^15),
#            mean(data.A[[p_score]]^15)))

#   # segments(10^mean(log10(data.I[[p_score]]^15)),          0, 10^mean(log10(data.I[[p_score]]^15)),          1)
#   # segments(10^mean(log10(data.I[[p_score]][inds_AU]^15)), 0, 10^mean(log10(data.I[[p_score]][inds_AU]^15)), 1)
#   # segments(10^mean(log10(data.A[[p_score]]^15)),          0, 10^mean(log10(data.A[[p_score]]^15)),          1)

#   points(c(10^mean(log10(data.I[[p_score]]^15)),
#            10^mean(log10(data.I[[p_score]][inds_AU]^15)),
#            10^mean(log10(data.A[[p_score]]^15))),
#         c(ecdf(data.I[[p_score]]^15)(10^mean(log10(data.I[[p_score]]^15))),
#           ecdf(data.I[[p_score]][inds_AU]^15)(10^mean(log10(data.I[[p_score]][inds_AU]^15))),
#           ecdf(data.A[[p_score]]^15)(10^mean(log10(data.A[[p_score]]^15)))),
#         col=c("black", "red", "red"),
#         pch = c(19, 20, 19))

#   axis(1, at=xs,
#     labels=FALSE,
# pos=0)

#     axis(1, at=xl,
#     labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
# pos=0)

#   axis(2, at=seq(0, 1, by=0.2), pos=10^-8, line = 1.5)



#   text(1e-2, 0.1, paste0("Effect: ",round(IndsEffect(condition)*100,1),"%"))
#   # text(1e-2, 0.2, paste0("Effect: lin average",round(pl_exp_by_AU*100,1),"%"))

#   title(xlab = "Per-nucleotide accessibility across from miRNA nt 1-15")
#   title(ylab = "Cumulative fraction")
#   legend(x = 1e-8, y = 1, cex = 0.9, legend = c("Input library",
#                     "AGO2-miR-1 library",
#                     "Input matched for AU context-score"),
#                     col = c("black", kEquilCondColors[as.character(condition)],
#                                      kEquilCondColors[as.character(condition)]),
#                     lty = c(1, 1, 2), bty="n")

#   effect.df <- data.frame(conditions = c(0.4, 1.26, 4, 12.6, 40),
#                           effect     = sapply(c(0.4, 1.26, 4, 12.6, 40),
#                                               function(condition) {
#                                                 IndsEffect(condition)
#                                                 }))
#   effect.df <<- effect.df
#   barplot(effect.df$effect,
#           names.arg=effect.df$conditions,
#           ylim=c(0, 1),
#           col=kEquilCondColors[as.character(effect.df$conditions)],
#           border=NA,
#           xlab="[AGO2-miR-1] (%v/v)")

#   }


# PlotSiteKdsVsRepression <- function(mirna, experiment, n_constant, sitelist,
#                                  cutoff=FALSE, bulk=FALSE, cat_colors=FALSE,
#                                  noncanon=FALSE, merge=FALSE, repression_df=repression.df,
#                                  xpos=20, ypos=20) {
#   dev.new(xpos=xpos, ypos=ypos , height=4, width=5)

#   par(kPlotParameters)

#   # mirna = "miR-1"
#   data <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
#   kds <- GetSiteKds(mirna, experiment, n_constant, sitelist)
#   data_seqs <- data[,1]
#   names(data_seqs) <- rownames(data)

#   site_type = "8mer"
#   sites_all <- unlist(unique(subset(repression_df,mir==mirna,select=site_type)))
#   sites_all <- sites_all[sites_all!="nosite"]
#   site_seqs <- sapply(sites_all, function(site) {
#     data_seqs[site]
#   })
#   names(site_seqs) <- sites_all
#   print(site_seqs)
#   site_seqs_noncanonical <- site_seqs
#   site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer"], site_seqs_noncanonical, invert=TRUE)]
#   site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer-A1"], site_seqs_noncanonical, invert=TRUE)]
#   site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer-m8"], site_seqs_noncanonical, invert=TRUE)]
#   site_seqs_noncanonical <- site_seqs_noncanonical[names(site_seqs_noncanonical)!="8mer-bG(6.7)"]
#   site_seqs_noncanonical <- site_seqs_noncanonical[names(site_seqs_noncanonical)!="7mer-m8bG(6.7)"]

#   site_seqs_canonical <- setdiff(sites_all,names(site_seqs_noncanonical))
#   print("canonical:")
#   print(site_seqs_canonical)
#   print("noncanonical:")
#   print(site_seqs_noncanonical)

#   # Option to merge all the noncanonical sites into a single site type:
#   if (merge == TRUE){
#     # Define the indeces of the noncanon
#     inds_keep <- which(repression_df$site_type %in% 
#                        names(site_seqs_noncanonical))
#     # Reassign the names of the noncanonical sites:
#     repression_df$site_type[inds_keep] <- "Noncanonical"
#     print(unique(repression_df$site_type))
#     sites_all <- c(site_seqs_canonical, "Noncanonical")
#   }
#   # Get the reduced matrix of fold change vaues:
#   out <- sapply(sites_all, function(site) {
#     print(site)
#     if (cutoff != FALSE &
#         (site %in% names(site_seqs_noncanonical) |
#           site == "Noncanonical")) {
#       reduced_frame <- subset(repression_df, 
#         mir==mirna & site_type==site & log_kd<=cutoff,
#         select=c(log_fc, log_kd))
#       if (site == "9mer-m13.21") {
#           print(reduced_frame)
#       }
#     } else {
#       reduced_frame <- subset(repression_df, 
#         mir == mirna & site_type == site,
#         select=c(log_fc, log_kd)
#       )      
#     }
#     print(dim(reduced_frame))
#     mean_values <- colMeans(reduced_frame)
#     sd_values <- apply(reduced_frame,2,sd)/sqrt(nrow(reduced_frame)-1)
#     print(sd_values)
#     if (bulk == TRUE) {
#       mean_values[2] <- log(kds[site,]$Mean,base=2)
#       sd_values[2] <- log(kds[site,]$Mean/(kds[site,]$Lower_CI),base=2)
#     }
#     return(c(mean_values,sd_values))
#   })

#   colnames(out) <- sites_all
#   print(out)
#   ind_remove <- !is.na(out[3,])
#   print(ind_remove)
#   out <- out[,ind_remove]
#   print(out)
#   if (noncanon == FALSE) {
#     out <- out[,c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8")]
#   }
#     sites_all <- colnames(out)


#   # print(out)
#       xmin <- 1
#     xmax <- 2000
#     ymin <- -1.4
#     ymax <- 0.4
#     nosite_rep <- mean(subset(repression_df,mir==mirna & site_type=="nosite",select=c(log_fc))[,1])
#   plot(1, type = "n",xlim=c(xmin, xmax*(xmax/xmin)^0.4), ylim = c(ymin,ymax),log='x',pch=19, lwd=2,ann=FALSE,axes=FALSE)
#   segments(xmin,0,xmax,0,lty=2)
#   arrows(1/(2^out[2,]), out[1,]+out[3,] - nosite_rep,1/(2^out[2,]),out[1,]-out[3,]-nosite_rep,length=0.05, lwd=1.5,angle=90, code=3)
#   arrows(1/(2^(out[2,]+out[4,])), out[1,]-nosite_rep,1/(2^(out[2,]-out[4,])),out[1,]-nosite_rep,length=0.05, lwd=1.5,angle=90, code=3)
#   if (cat_colors == TRUE){
#     colors.sites = kSiteCategoryColors[sites_all,]
#   } else {
#     colors.sites = kSiteColors[sites_all,]
#   }

#   points(1/(2^out[2,]),out[1,]-nosite_rep,col="black", bg=colors.sites,pch=21,cex=1.5)


#   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   ys <- seq(ymin,ymax,by=0.1)

#   xs <- xs[xs >= xmin & xs <= xmax]
#   ys <- ys[ys >= ymin & ys <= ymax]

#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   yl <- seq(ymin,ymax, by= 0.2)

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin, lwd=2)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=round(yl,2),
#        pos=xmin, las=2, lwd=0, )
#   axis(2, at=ys, labels=FALSE,
#        pos=xmin, lwd=2)


#   remove <- unique(which(is.na(out[1,])), which(is.na(out[1,])))
#   if (length(remove) > 0) {
#     out <- out[,-remove]  
#   }
#   text(x=800, y=0.3, mirna,cex=0.9)

#   text(x=800, y=0.2, eval(substitute(expression(italic(r) == x), 
#             list(x = round(cor(out[2,],out[1,]),3)))),cex=0.9)
#   if (cutoff != FALSE) {
#   text(x=500, y=0.1, eval(substitute(expression(log[2](italic(K)[D][italic(nc)]) <= x),
#             list(x = cutoff))), cex = 0.9)
#   }
#   sites_all <- sites_all[order(out[2,])]
#   if (cat_colors == TRUE){
#     colors.sites = kSiteCategoryColors[sites_all,]
#   } else {
#     colors.sites = kSiteColors[sites_all,]
#   }
#   legend(x=2000,y=0.4, legend=sites_all, cex=0.8,bty="n", pch=19,col=colors.sites, ncol=1)
#   sites_all <<- sites_all
#   title(xlab=expression(italic(K)[D]))
#   title(ylab=expression(log[2](paste("fold change"))))
# }


# PlotAllSiteKdsVsRepression <- function(experiment, n_constant, sitelist,
#                                  cutoff=FALSE, bulk=FALSE,
#                                  noncanon=FALSE, repression_df = repression.df, merge=FALSE, xpos=20, 
#                                  ypos=20) {
#   dev.new(xpos=xpos, ypos=ypos , height=5, width=5)

#   par(kPlotParameters)
#     xmin <- 1/2000
#     xmax <- 1
#     ymin <- -1.4
#     ymax <- 0.4
#   plot(1, type = "n",xlim=c(xmax, xmin), ylim = c(ymin,ymax),log='x',pch=19, lwd=2,ann=FALSE,axes=FALSE)

#   for (mirna in c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")) {
#   data <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
#   kds <- GetSiteKds(mirna, experiment, n_constant, sitelist)
#   print(kds)
#   data_seqs <- data[,1]
#   names(data_seqs) <- rownames(data)

#   print(data_seqs)
#   site_type = "8mer"
#   sites_all <- unlist(unique(subset(repression_df,mir==mirna,select=site_type)))
#   print(sites_all)
#   sites_all <- sites_all[sites_all!="nosite"]
#   site_seqs <- sapply(sites_all, function(site) {
#     data_seqs[site]
#   })
#   names(site_seqs) <- sites_all
#   print('sites_all')
#   print(sites_all)

#   print(site_seqs)
#   print(data_seqs)
#   site_seqs_noncanonical <- site_seqs
#   site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer"], site_seqs_noncanonical, invert=TRUE)]
#   site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer-A1"], site_seqs_noncanonical, invert=TRUE)]
#   site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer-m8"], site_seqs_noncanonical, invert=TRUE)]
#   site_seqs_noncanonical <- site_seqs_noncanonical[names(site_seqs_noncanonical)!="8mer-bG(6.7)"]
#   site_seqs_noncanonical <- site_seqs_noncanonical[names(site_seqs_noncanonical)!="7mer-m8bG(6.7)"]

#   site_seqs_canonical <- setdiff(sites_all,names(site_seqs_noncanonical))
#   print("canonical")
#   print(site_seqs_canonical)
#   print("noncanonical")
#   print(site_seqs_noncanonical)
#   if (merge == TRUE){
#     repression_df$site_type[which(repression_df$site_type %in% names(site_seqs_noncanonical))] <- "Noncanonical"
#       print(unique(repression_df$site_type))
#     sites_all <- c(site_seqs_canonical, "Noncanonical")
#   }

#   out <- sapply(sites_all, function(site) {
#     print(site)
#     if (cutoff != FALSE &
#         (site %in% names(site_seqs_noncanonical) |
#           site == "Noncanonical")) {
#       reduced_frame <- subset(repression_df, 
#         mir==mirna & site_type==site & log_kd<=cutoff,
#         select=c(log_fc, log_kd))
#       if (site == "9mer-m13.21") {
#           print(reduced_frame)
#       }
#     } else {
#       reduced_frame <- subset(repression_df, 
#         mir == mirna & site_type == site,
#         select=c(log_fc, log_kd)
#       )      
#     }
#     print(dim(reduced_frame))
#     mean_values <- colMeans(reduced_frame)
#     sd_values <- apply(reduced_frame,2,sd)/sqrt(nrow(reduced_frame)-1)
#     print(sd_values)
#     if (bulk == TRUE) {
#       mean_values[2] <- log(kds[site,]$Mean,base=2)
#       sd_values[2] <- log(kds[site,]$Mean/(kds[site,]$Lower_CI),base=2)
#     }
#     return(c(mean_values,sd_values))
#   })

#   colnames(out) <- sites_all
#   print(out)
#   ind_remove <- !is.na(out[3,])
#   print(ind_remove)
#   out <- out[,ind_remove]
#   print(out)
#   if (noncanon == FALSE) {
#     out <- out[,c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8")]
#     sites_all <- colnames(out)
#   }
#   nosite_rep <- mean(subset(repression_df,
#                             mir==mirna & site_type=="nosite",
#                             select=c(log_fc))[, 1])
#   if (mirna == "let-7a") {
#     color <- "gray70"
#     kMirnaColors[mirna] <- "gray70"
#   } else {
#     color <- kMirnaColors[mirna]    
#   }
#   arrows(2^out[2,],
#          out[1,]+out[3,] - nosite_rep,
#          2^out[2,],
#          out[1,] - out[3,] - nosite_rep, length=0.05, lwd=1, angle=90,
#          col=color, code=3)
#   arrows(2^(out[2,] + out[4,]),
#          out[1,] - nosite_rep,
#          2^(out[2,] - out[4,]),
#          out[1,] - nosite_rep, length=0.05, lwd=1, angle=90, col=color,
#          code=3)
#   points(2^out[2,],
#          out[1,] - nosite_rep,
#          col=color,
#          pch=19)

#   remove <- unique(which(is.na(out[1,])), which(is.na(out[1,])))
#   if (length(remove) > 0) {
#     out <- out[,-remove]  
#   }
#   sites_all <- sites_all[order(out[2,])]
#   if (mirna == "let-7a") {
#     all.sites <- out
#   } else if (mirna == "miR-1") {
#     exclude.l7.sites <- out

#   } else {
#     all.sites <- cbind(all.sites, out)
#     exclude.l7.sites <- cbind(exclude.l7.sites, out)
#   }
#   }

#   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   print(xs)
#   ys <- seq(ymin,ymax,by=0.1)

#   xs <- xs[xs >= xmin & xs <= xmax]
#   ys <- ys[ys >= ymin & ys <= ymax]

#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   yl <- seq(ymin,ymax, by= 0.2)

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin, lwd=2)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=round(yl,2),
#        pos=xmax, las=2, lwd=0, )
#   axis(2, at=ys, labels=FALSE,
#        pos=xmax, lwd=2)

#   legend(x=0.8,y=-0.8, legend=mirnas.all, bty="n", pch=19,col=kMirnaColors[mirnas.all], ncol=1)
#   sites_all <<- sites_all

#   text(x=1e-3, y=0, eval(substitute(expression(italic(r) == x), 
#             list(x = round(cor(all.sites[2,],all.sites[1,]),3)))),
#             col="gray50")
#   text(x=1e-3, y=-0.07, eval(substitute(expression(italic(r) == x), 
#             list(x = round(cor(exclude.l7.sites[2,],exclude.l7.sites[1,]),3)))))
#   if (cutoff != FALSE) {
#   text(x=1/500, y=0.1, eval(substitute(expression(log[2](italic(K)[D][italic(nc)]) <= x),
#             list(x = cutoff))), cex = 0.9)
#   }

#   fit <- lm(all.sites[1, ] ~ log10(2^all.sites[2, ]))
#   m <- fit$coefficients[2]
#   b <- fit$coefficients[1]
#   x_line <- 10^(seq(log10(xmin), log10(xmax),length = 20))
#   y_line <- m*log10(x_line) + b
#   lines(x_line, y_line, lty = 2,lwd = 1,col="gray50")
  
#   fit <- lm(exclude.l7.sites[1, ] ~ log10(2^exclude.l7.sites[2, ]))
#   m <- fit$coefficients[2]
#   b <- fit$coefficients[1]
#   x_line <- 10^(seq(log10(xmin), log10(xmax),length = 20))
#   y_line <- m*log10(x_line) + b
#   lines(x_line, y_line, lty = 2,lwd = 1)
#   # lowess_fit <<- lowess(log_fc ~ log_kd,
#   #                      data=subset(repression.df,
#   #                                  site_type %in% canonical.sites,
#   #                                  select = c(log_kd,log_fc)))
#   # print(lowess_fit)
#   # lines(2^lowess_fit$x,lowess_fit$y)
#   # lowess_fit <<- lowess(log_fc ~ log_kd,
#   #                      data=subset(repression.df,
#   #                                  mir != "let-7a" & site_type %in% canonical.sites,
#   #                                  select = c(log_kd,log_fc)))
#   # print(lowess_fit)
#   # lines(2^lowess_fit$x,lowess_fit$y)

#   # lowess_fit <<- lowess(log_fc ~ log_kd,
#   #                      data=subset(repression.df,
#   #                                  mir != "let-7a" & site_type == "8mer",
#   #                                  select = c(log_kd,log_fc)))
#   # print(lowess_fit)
#   # lines(2^lowess_fit$x,lowess_fit$y)


#   if (merge == TRUE) {
#   text(x=1/500, y=-0.09, eval(substitute(expression(log[2](italic(K)[D][italic(noncanon)]) <= x),
#             list(x = cutoff))))
#   }


#   title(xlab=expression(italic(K)[D]))
#   title(ylab=expression(log[2](paste("fold change"))))
#   # text(x=30,y=0.1,round(cor(out[2,],out[1,]),3),cex=1.5)

# }

# CompareTwoSiteTypesForRepression <- function(n_constant, sitelist, sample_text,
#                                  num_boxes=10,
#                                  cutoff=FALSE, bulk=FALSE,
#                                  repression_df = repression.df, merge=FALSE, xpos=20, 
#                                  ypos=20) {
#   dev.new(xpos=xpos, ypos=ypos , height=6, width=10)
#   par(kPlotParameters)
#   canonical.sites <- c("8mer", "7mer-m8", "7mer-A1", "6mer")
#   par(mfrow=c(2, 3))
#   for (mirna in mirnas.all) {
#     rep.df <- subset(repression_df, mir==mirna & site_type %in% canonical.sites)
#     round_kd <- round(rep.df$log_kd)    
#     rep.df$breaks <- round_kd
#     print(breaks)
#     breaks <<- breaks
#     print(unique(breaks))
#     rep.df$breaksbysite <- paste0(rep.df$site_type,"|", rep.df$breaks)
#     print(unique(rep.df$breaksbysite))
#     new_matrix <- aggregate(log_fc ~ site_type + breaks, data=rep.df, mean)
#     print(new_matrix)
#     output <- matrix(0,nrow=length(unique(rep.df$site_type)), ncol = length(unique(rep.df$breaks)))
#     rownames(output) <- canonical.sites
#     print(output)
#     print("hi")
#     print(sort(unique(breaks)))
#     colnames(output) <- sort(unique(round_kd))
#     print("hi")

#     print(output)
#     for (row in seq(nrow(new_matrix))) {
#       print(row)
#       print(output)
#       print(new_matrix[row,])
#       output[as.character(new_matrix[row,1]), as.character(new_matrix[row, 2])] <- new_matrix[row, 3]
#     }
#     print(output)
#     barplot(output, beside=TRUE,col=kSiteColors[canonical.sites,],ylim = c(-2, 0.5))
#     text(5, 0.2, mirna)
# }

#     text(10, -1.9, sample_text)
# }
# # CompareTwoSiteTypesForRepression(5, "paper", "new repression", num_boxes=check, xpos = 20, ypos = 20)
# # dev.copy2pdf(file = "2017_Paper/Figure4G1_raw_v8.pdf")

# # repression.df <- repression_old.df
# # CompareTwoSiteTypesForRepression(5, "paper", "old repression", num_boxes=check, xpos = 20, ypos = 20)
# # dev.copy2pdf(file = "2017_Paper/Figure4G2_raw_v8.pdf")




# PlotSingleFlankKdsVsRepression <- function(mirna, experiment, n_constant,
#                                            sitelist, site,
#                                            cutoff=FALSE, bulk=FALSE,
#                                            noncanon=FALSE, repression_df = repression.df, merge=FALSE, xpos=20, 
#                                            ypos=20) {
#   dev.new(xpos=xpos, ypos=ypos , height=5, width=5)
#   par(kPlotParameters)
#     xmin <- 1/10000
#     xmax <- 10
#     ymin <- -4
#     ymax <- 2
#   plot(1, type = "n", 
#        xlim=c(xmax, xmin),
#        ylim=c(ymin, ymax),
#        log='x',
#        pch=19,
#        lwd=2,
#        ann=FALSE,
#        axes=FALSE)
#   segments(xmin, 0, xmax, 0, lty = 2)
#   # Get the flanking kds for the site:
#   # kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
#   # print(kds)
#   # Subset the matrix for only the mirna and site type:
#     if (mirna == "all"){
#   rep.df <- data.frame(site_sequence = c(), log_fc = c(), log_kd = c(), site_type = c())
#    for (mirna in mirnas.all) {
#            rep.df_temp <- subset(repression_df,
#     mir == mirna & site_type == site, select=c(site_sequence, log_fc, log_kd, site_type))   
#            rep.df_temp$log_fc <- rep.df_temp$log_fc - mean(unlist(subset(repression_df, mir==mirna & site == "nosite",select=c(log_fc))))
#     rep.df <- rbind(rep.df, rep.df_temp)
#    }
#    } else {
#       rep.df <- subset(repression_df,
#     mir==mirna & site_type == site, select=c(site_sequence, log_fc, log_kd, site_type))

#    }



#   # Convert the sequence on the left into the flanking dinucleotides,
#   # for the coloration of the points.
#   ConvertSequenceToFlanks <- function(sequence){
#     left <- substr(sequence, 1, 2)
#     right <- substr(sequence, nchar(sequence) - 1, nchar(sequence))
#     out <- paste0(left, right, collapse="")
#     return(out)
#   }
#   rep.df$flanks <- sapply(rep.df$site_sequence, ConvertSequenceToFlanks)
#   agg_fc <- aggregate(log_fc ~ flanks, rep.df, function(x) {
#     c(mean = mean(c(x)),
#       se = sd(c(x))/sqrt(length(c(x))-1))
#     })
#   agg_kd <- aggregate(log_kd ~ flanks, rep.df, function(x) {
#     c(mean = mean(x))
#     })
#   agg_kd <<- agg_kd
#   rep.df <<- rep.df
#   agg_fc <<- agg_fc
#   # Make the color vector:
#   colors <- sapply(rep.df$flanks, GetColorFunction, alpha=0.9)
#   colors_sum <- sapply(agg_fc$flanks, GetColorFunction)
#   points(2^rep.df$log_kd, rep.df$log_fc, pch=19, col = colors)
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   ys <- seq(ymin,ymax,by=0.1)

#   xs <- xs[xs >= xmin & xs <= xmax]
#   ys <- ys[ys >= ymin & ys <= ymax]

#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   yl <- seq(ymin,ymax, by= 0.5)

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin, lwd=2)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=round(yl,2),
#        pos=xmax, las=2, lwd=0, )
#   axis(2, at=ys, labels=FALSE,
#        pos=xmax, lwd=2)
#   title(xlab=expression(italic(K)[D]))
#   title(ylab=expression(log[2](paste("fold change"))))
#   text(x = 10^-3, y = 1.5, mirna)
#   text(x = 10^-3, y = 1.2, site)

#   # plot(1, type = "n", 
#   #      xlim=c(xmax, xmin),
#   #      ylim=c(ymin, ymax),
#   #      log='x',
#   #      pch=19,
#   #      lwd=2,
#   #      ann=FALSE,
#   #      axes=FALSE)
#   # segments(xmin, 0, xmax, 0, lty = 2)

#   # points(2^agg_kd[,2], agg_fc$log_fc[,1], col = colors_sum)


#   # axis(1, at=xl,
#   #      labels=sapply(xl, function(name) {
#   #        eval(substitute(expression(10^x), list(x=log10(name))))
#   #      }),
#   #      pos=ymin, lwd=0)
#   # axis(1, at=xs, labels=FALSE,
#   #      pos=ymin, lwd=2)
#   # # Label the axis at each order of magnitude.

#   # axis(2, at=yl,
#   #      labels=round(yl,2),
#   #      pos=xmax, las=2, lwd=0, )
#   # axis(2, at=ys, labels=FALSE,
#   #      pos=xmax, lwd=2)


# }


# PlotAllCanonicalKdsVsRepression <- function(mirna, experiment, n_constant,
#                                            sitelist, 
#                                            cutoff=FALSE, bulk=FALSE,
#                                            noncanon=FALSE, repression_df = repression.df, merge=FALSE, xpos=20, 
#                                            ypos=20) {
#   dev.new(xpos=xpos, ypos=ypos , height=5, width=5)
#   par(kPlotParameters)
#     xmin <- 1/10000
#     xmax <- 10
#     ymin <- -4
#     ymax <- 2
#   plot(1, type = "n", 
#        xlim=c(xmax, xmin),
#        ylim=c(ymin, ymax),
#        log='x',
#        pch=19,
#        lwd=2,
#        ann=FALSE,
#        axes=FALSE)
#   segments(xmin, 0, xmax, 0, lty = 2)
#   # Get the flanking kds for the site:
#   # Subset the matrix for only the mirna and site type:
#   if (mirna == "all"){
#    rep.df <- subset(repression_df,
#     site_type %in% canonical.sites, select=c(site_sequence, log_fc, log_kd, site_type))
   
#    } else {
#       rep.df <- subset(repression_df,
#     mir==mirna & site_type %in% canonical.sites, select=c(site_sequence, log_fc, log_kd, site_type))

#    }
#   # Convert the sequence on the left into the flanking dinucleotides,
#   # for the coloration of the points.
#   ConvertSequenceToFlanks <- function(sequence){
#     left <- substr(sequence, 1, 2)
#     right <- substr(sequence, nchar(sequence) - 1, nchar(sequence))
#     out <- paste0(left, right, collapse="")
#     return(out)
#   }
#   rep.df$flanks <- sapply(rep.df$site_sequence, ConvertSequenceToFlanks)
#   agg_fc <- aggregate(log_fc ~ flanks, rep.df, function(x) {
#     c(mean = mean(c(x)),
#       se = sd(c(x))/sqrt(length(c(x))-1))
#     })
#   agg_kd <- aggregate(log_kd ~ flanks, rep.df, function(x) {
#     c(mean = mean(x))
#     })
#   agg_kd <<- agg_kd
#   rep.df <<- rep.df
#   agg_fc <<- agg_fc
#   # Make the color vector:
#   colors.site <- sapply(kSiteColors[rep.df$site_type,],function(name){
#     rgb_values <- c(col2rgb(name))/255
#     red <- rgb_values[1]
#     green <- rgb_values[2]
#     blue <- rgb_values[3]

#     return(rgb(red, green, blue ,alpha=0.2))})

#     colors.flanks <- sapply(rep.df$flanks, GetColorFunction, alpha=0.1)

#   colors_sum <- sapply(agg_fc$flanks, GetColorFunction)
#   points(2^rep.df$log_kd, rep.df$log_fc, pch=19, col = colors.site)
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   print(xs)
#   ys <- seq(ymin,ymax,by=0.1)

#   xs <- xs[xs >= xmin & xs <= xmax]
#   ys <- ys[ys >= ymin & ys <= ymax]

#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   yl <- seq(ymin,ymax, by= 0.5)

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin, lwd=2)
#   # Label the axis at each order of magnitude.
#   text(x=1, y=-3.5, eval(substitute(expression(italic(r) == x), 
#             list(x = round(cor(rep.df$log_kd,rep.df$log_fc),3)))))

#   axis(2, at=yl,
#        labels=round(yl,2),
#        pos=xmax, las=2, lwd=0, )
#   axis(2, at=ys, labels=FALSE,
#        pos=xmax, lwd=2)
#   title(xlab=expression(italic(K)[D]))
#   title(ylab=expression(log[2](paste("fold change"))))
#   text(x = 10^-3, y = 1.5, mirna)
#   for (site in canonical.sites) {
#     lowess_object <- lowess(log_fc ~ log_kd, subset(rep.df, site_type==site),f=1)
#     lines(2^lowess_object$x,lowess_object$y, col=kSiteColors[site,])
#   }
#   # plot(1, type = "n", 
#   #      xlim=c(xmax, xmin),
#   #      ylim=c(ymin, ymax),
#   #      log='x',
#   #      pch=19,
#   #      lwd=2,
#   #      ann=FALSE,
#   #      axes=FALSE)
#   # segments(xmin, 0, xmax, 0, lty = 2)


#   # points(2^rep.df$log_kd, rep.df$log_fc, pch=1, col = colors.flanks)
#   #   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   # print(xs)
#   # ys <- seq(ymin,ymax,by=0.1)

#   # xs <- xs[xs >= xmin & xs <= xmax]
#   # ys <- ys[ys >= ymin & ys <= ymax]

#   # xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   # yl <- seq(ymin,ymax, by= 0.5)

#   # axis(1, at=xl,
#   #      labels=sapply(xl, function(name) {
#   #        eval(substitute(expression(10^x), list(x=log10(name))))
#   #      }),
#   #      pos=ymin, lwd=0)
#   # axis(1, at=xs, labels=FALSE,
#   #      pos=ymin, lwd=2)
#   # # Label the axis at each order of magnitude.

#   # axis(2, at=yl,
#   #      labels=round(yl,2),
#   #      pos=xmax, las=2, lwd=0, )
#   # axis(2, at=ys, labels=FALSE,
#   #      pos=xmax, lwd=2)
#   #   xmin <- 1/2000
#   #   xmax <- 10
#   #   ymin <- -1
#   #   ymax <- 0.5

#   # plot(1, type = "n", 
#   #      xlim=c(xmax, xmin),
#   #      ylim=c(ymin, ymax),
#   #      log='x',
#   #      pch=19,
#   #      lwd=2,
#   #      ann=FALSE,
#   #      axes=FALSE)
#   # segments(xmin, 0, xmax, 0, lty = 2)

#   # points(2^agg_kd[,2], agg_fc$log_fc[,1], col = colors_sum)


#   # axis(1, at=xl,
#   #      labels=sapply(xl, function(name) {
#   #        eval(substitute(expression(10^x), list(x=log10(name))))
#   #      }),
#   #      pos=ymin, lwd=0)
#   # axis(1, at=xs, labels=FALSE,
#   #      pos=ymin, lwd=2)
#   # # Label the axis at each order of magnitude.

#   # axis(2, at=yl,
#   #      labels=round(yl,2),
#   #      pos=xmax, las=2, lwd=0, )
#   # axis(2, at=ys, labels=FALSE,
#   #      pos=xmax, lwd=2)


# }

# PlotLinearModelKdVsRepressionCanonical <- function(
#   experiment, n_constant, sitelist, cutoff=FALSE, bulk=FALSE, noncanon=FALSE,
#   merge=FALSE, repression_df=repression.df, xpos=20, ypos=20) {
#   # Get the flanking kds for the site:
#   # Subset the matrix for only the mirna and site type:
#    rep.df <- subset(repression_df,
#                     site_type %in% canonical.sites,
#                     select=c(mir, site_sequence, log_fc, log_kd, site_type))
#   # Convert the sequence on the left into the flanking dinucleotides,
#   # for the coloration of the points.
#   # Define subfunction:
#   ConvertSequenceToFlanks <- function(sequence) {
#     left  <- substr(sequence,
#                     1,
#                     2)
#     right <- substr(sequence,
#                     nchar(sequence) - 1,
#                     nchar(sequence))
#     out   <- paste0(left,
#                     right,
#                     collapse="")
#     return(out)
#   }
#   # Use function.
#   rep.df$flank <- sapply(rep.df$site_sequence,
#                          ConvertSequenceToFlanks)
#   # Solve all four linear models:
#   attach(rep.df)
#   lm.kd  <- lm(log_kd ~ mir + site_type + flank)
#   lm.kd2 <- lm(log_kd ~ mir * site_type + flank)
#   lm.fc  <- lm(log_fc ~ mir + site_type + flank)
#   lm.fc2 <- lm(log_fc ~ mir * site_type + flank)
#   detach(rep.df)
#   # Get colors for plot:
#   colors.site <- sapply(kSiteColors[rep.df$site_type,],
#                         ConvertRColortoRGB,
#                         alpha = 0.1)
#   colors.flank <- sapply(rep.df$flank,
#                          GetColorFunction,
#                          alpha=0.1)
#   colors.mat <- cbind(colors.site, colors.flank)
#   colors.site <<- colors.site
#   colors.flank <<- colors.flank
#   colors.mat <<- colors.mat
#   colors_sum <- sapply(agg_fc$flank, GetColorFunction)
#   # points(2^model.kds,2^rep.df$log_kd, pch=1, col = colors.site)

#   dev.new(xpos   = xpos,
#           ypos   = ypos,
#           height = 8,
#           width  = 12)
#   par(kPlotParameters)
#   par(mfrow=c(3, 4))
#   xmin <- 1/10000
#   xmax <- 10
#   ymin <- -4.5
#   ymax <- 1.5
#   for (lm.ind in seq(2)) {
#     for (color.ind in seq(2)) {
#     x_data <- fitted(list(lm.kd, lm.kd2)[[lm.ind]])
#     y_data <- rep.df$log_kd
#     plot(2^x_data,
#          2^y_data,
#          xlim = c(xmax, xmin),
#          ylim = c(xmax, xmin),
#          log  = 'xy',
#          pch  = 19,
#          col  = colors.mat[,color.ind],
#          lwd  = 2,
#          ann  = FALSE,
#          axes = FALSE)
#     segments(xmin,
#              0,
#              xmax,
#              0,
#              lty = 2)
#     xs <- c(sapply(seq(floor(log10(xmin)),
#                        ceiling(log10(xmax))),
#                    function(x) seq(10)*10^x))
#     ys <- seq(ymin,
#               ymax,
#               by=0.1)
#     xs <- xs[xs >= xmin & xs <= xmax]
#     ys <- ys[ys >= ymin & ys <= ymax]

#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#     yl <- seq(ymin,ymax, by= 0.5)

#     axis(1,
#          at     = xl,
#          labels = sapply(xl,
#                          function(name) {
#            eval(substitute(expression(10^x),
#                            list(x=log10(name))))
#            }),
#          pos    =xmax,
#          lwd    = 0)
#     axis(1,
#          at     = xs,
#          labels = FALSE,
#          pos    = xmax,
#          lwd    = 2)
#     # Label the axis at each order of magnitude.
#     text(x=1e-3,
#          y=0.5,
#          eval(substitute(expression(italic(r^2) == x), 
#               list(x = round(cor(x_data, y_data)^2,3)))))
#     axis(2,
#          at     = xl,
#          labels = sapply(xl,
#                          function(name) {
#            eval(substitute(expression(10^x),
#                            list(x=log10(name))))
#            }),
#          pos    =xmax,
#          lwd    = 0)
#     axis(2,
#          at     = xs,
#          labels = FALSE,
#          pos    = xmax,
#          lwd    = 2)
#   title(xlab=expression(italic(K)[D],predicted))
#   title(ylab=expression(italic(K)[D],measured))

# }}

#   sapply(seq(2), function(lm.ind) {

#   sapply(seq(2), function(color.ind) {
#     x_data <- fitted(list(lm.fc, lm.fc2)[[lm.ind]])
#     y_data <- rep.df$log_fc
#     plot(0, type = "n",
#          xlim = c(ymax, ymin),
#          ylim = c(ymax, ymin),
#          pch  = 19,
#          ann  = FALSE,
#          axes = FALSE)
#     segments(ymin,
#              0,
#              ymax,
#              0,
#              lty = 2)
#     segments(0,
#              ymin,
#              0,
#              ymax,
#              lty = 2)
#     points(x_data,
#            y_data,
#            col = colors.mat[,color.ind],
#            pch=19)

#     ys <- seq(ymin,
#               ymax,
#               by=0.5)
#     ys <- ys[ys >= ymin & ys <= ymax]

#     yl <- seq(ymin, ymax, by= 1)

#     axis(1,
#          at     = yl,
#          pos    =ymax,
#          lwd    = 0)
#     axis(1,
#          at     = ys,
#          labels = FALSE,
#          pos    = ymax,
#          lwd    = 2)
#     # Label the axis at each order of magnitude.
#     text(x=-3,
#          y=0.5,
#          eval(substitute(expression(italic(r^2) == x), 
#               list(x = round(cor(x_data, y_data)^2,3)))))
#     axis(2,
#          at     = yl,
#          pos    =ymax,
#          lwd    = 0)
#     axis(2,
#          at     = ys,
#          labels = FALSE,
#          pos    = ymax,
#          lwd    = 2)
#   title(xlab=expression(log[2](paste("fold change"))))

#   title(ylab=expression(log[2](paste("fold change"))))

# })})


#     GetPredictedValuesSiteMiR <- function(model,error=FALSE){
#       if (error == TRUE) {
#         extract <- "se.fit"
#       } else {
#         extract <- "fit"
#       }
#       sapply(canonical.sites, function(site) {
#         sapply(mirnas.all, function(mirna) {
#           prediction <- predict(model,data.frame(mir=c(mirna), site_type=c(site),flank=c("AAAA")),se.fit=TRUE)
#           # predict(lm.kd,data.frame(mir=c("miR-1"), site_type=c("8mer"),flank=c("AAAA")),se.fit=TRUE)
#           return(c(prediction[[extract]]))
#           })
#       })
#     }

#     GetPredictedValuesFlanks <- function(model,error=FALSE){
#       if (error == TRUE) {
#         extract <- "se.fit"
#       } else {
#         extract <- "fit"
#       }
#       sapply(sort(unique(rep.df$flank)), function(flank) {
#           prediction <- predict(model,data.frame(mir=c("let-7a"), site_type=c("6mer"),flank=c(flank)),se.fit=TRUE)
#           # predict(lm.kd,data.frame(mir=c("miR-1"), site_type=c("8mer"),flank=c("AAAA")),se.fit=TRUE)
#           return(c(prediction[[extract]]))
#           })
#     }


#     y.ms.kd <<- GetPredictedValuesSiteMiR(lm.kd2)
#     e.ms.kd <<- GetPredictedValuesSiteMiR(lm.kd2, error=TRUE)
#     y.ms.fc <<- GetPredictedValuesSiteMiR(lm.fc2)
#     e.ms.fc <<- GetPredictedValuesSiteMiR(lm.fc2, error=TRUE)

#     y.f.kd <<- GetPredictedValuesFlanks(lm.kd2)
#     e.f.kd <<- GetPredictedValuesFlanks(lm.kd2, error=TRUE)
#     y.f.fc <<- GetPredictedValuesFlanks(lm.fc2)
#     e.f.fc <<- GetPredictedValuesFlanks(lm.fc2, error=TRUE)




#     color.flanks <- sapply(sort(unique(rep.df$flank)),GetColorFunction,alpha=0.5)
#     color.site <- kSiteColors[rep(canonical.sites,each=5),]
#     plot(1, type="n", xlim=c(-12, -2), ylim=c(-1.5, 1))
#     arrows(y.ms.kd,
#       y.ms.fc-e.ms.fc,
#       y.ms.kd,
#       y.ms.fc+e.ms.fc,
# ,length=0.07, angle=90, code=3,lwd=0.5)

#     arrows(y.ms.kd - e.ms.kd,
#       y.ms.fc,
#       y.ms.kd + e.ms.kd,
#       y.ms.fc,
# ,length=0.07, angle=90, code=3,lwd=0.5)


#     kd_rep_mirsite <<- lm(c(y.ms.fc) ~ c(y.ms.kd), weights=1/c(e.ms.fc))
#     kd_rep_flank <<- lm(c(y.f.fc) ~ c(y.f.kd), weights=1/c(e.f.fc))

#     abline(kd_rep_flank,lty=2,col="gray")
#     abline(kd_rep_mirsite,lty=2)

#     points(y.ms.kd, y.ms.fc,col=color.site,pch=19)

#     axis(1,
#          at     = xl,
#          labels = sapply(xl,
#                          function(name) {
#            eval(substitute(expression(10^x),
#                            list(x=log10(name))))
#            }),
#          pos    =xmax,
#          lwd    = 0)
#     axis(1,
#          at     = xs,
#          labels = FALSE,
#          pos    = xmax,
#          lwd    = 2)

#     axis(2,
#          at     = yl,
#          pos    =ymax,
#          lwd    = 0)
#     axis(2,
#          at     = ys,
#          labels = FALSE,
#          pos    = ymax,
#          lwd    = 2)
#   title(xlab=expression(italic(K)[D],predicted))
#   title(ylab=expression(log[2](paste("fold change"))))



#     plot(1, type="n", xlim=c(-12, -2), ylim=c(-1.5, 1))

#     arrows(y.f.kd,
#       y.f.fc-e.f.fc,
#       y.f.kd,
#       y.f.fc+e.f.fc,
# ,length=0.07, angle=90, code=3,lwd=0.2)

#     arrows(y.f.kd - e.f.kd,
#       y.f.fc,
#       y.f.kd + e.f.kd,
#       y.f.fc,
# ,length=0.07, angle=90, code=3,lwd=0.2)

#       abline(kd_rep_flank,lty=2)
#       abline(kd_rep_mirsite,lty=2,col="gray")

#     points(y.f.kd, y.f.fc,col=color.flanks,pch=19)
#   title(xlab=expression(italic(K)[D],predicted))

#   title(ylab=expression(log[2](paste("fold change"))))

#     print(kd_rep_mirsite)
#     print(kd_rep_flank)
# }


# # FOR PAPER
# PlotpairwiseKds <- function(mirna, experiment, start, stop, sitelist, combined.input=TRUE) {
#   kds.1 <- GetKds(mirna, experiment, start, stop, sitelist, combined.input = combined.input)
#   kds.2 <- GetKds(mirna, experiment, start, stop, sitelist, log.residual = TRUE, combined.input=combined.input)
#   par(kPlotParameters)
#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2))))))
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#     plot(kds.1,
#       kds.2,
#       col = kSiteColors[c(names(kds.1), "bg", "Ago"),],
#       pch = c(rep(20, length(kds.1) - 2), 1, 1),
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#     segments(xmin, xmin, xmax, xmax, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = paste0("Multinomial optimized parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0("Log-normal optimized parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

# }

# PlotpairwiseKdsInput <- function(mirna, experiment, start, stop, sitelist, log.residual=FALSE) {
#   kds.1 <- GetKds(mirna, experiment, start, stop, sitelist, log.residual=log.residual)
#   kds.2 <- GetKds(mirna, experiment, start, stop, sitelist, combined.input = FALSE, log.residual=log.residual)
#   par(kPlotParameters)
#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2))))))
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(kds.1,
#       kds.2,
#       col = kSiteColors[c(names(kds.1), "bg", "Ago"),],
#       pch = c(rep(20, length(kds.1) - 2), 1, 1),
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#     segments(xmin, xmin, xmax, xmax, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = paste0("Combined input parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0("Single input parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

# }

# PlotpairwiseFlankKds <- function(mirna, experiment, start, stop, site,
#                                  sitelist, combined.input=TRUE) {
#   kds.1 <- GetFlankKds(mirna, experiment, start, stop, site, sitelist,
#                        combined.input=combined.input)
#   kds.2 <- GetFlankKds(mirna, experiment, start, stop, site, sitelist,
#                        combined.input=combined.input, log.residual=TRUE)
#   par(kPlotParameters)
#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2))))))
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(kds.1,
#       kds.2,
#       col = sapply(names(kds.1), GetColorFunction),
#       pch = 20,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#     segments(xmin, xmin, xmax, xmax, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = paste0("Multinomial optimized parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0("Log-normal optimized parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

# }


# PlotpairwiseFlankKdsInput <- function(mirna, experiment, start, stop, site,
#                                       sitelist, log.residual=FALSE) {
#   print("622")
#   kds.1 <- GetFlankKds(mirna, experiment, start, stop, site, sitelist, log.residual=log.residual)
#   print("624")
#   kds.2 <- GetFlankKds(mirna, experiment, start, stop, site, sitelist, combined.input = FALSE, log.residual=log.residual)
#   print(kds.1)
#   print(kds.2)
#   par(kPlotParameters)
#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2))))))
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(kds.1,
#       kds.2,
#       col = sapply(names(kds.1), GetColorFunction),
#       pch = 20,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#     segments(xmin, xmin, xmax, xmax, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = paste0("Combined input parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0("Single input parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

# }





# PlotpairwiseFlankKdsSites <- function(mirna, experiment, start, stop, site1, site2, sitelist, combined.input=TRUE, log.residual=FALSE,position=FALSE,nosite=TRUE) {
#   kds.1 <- GetFlankKds(mirna, experiment, site1, start, stop, sitelist, combined.input = combined.input, log.residual=log.residual)
#   kds.2 <- GetFlankKds(mirna, experiment, site2, start, stop, sitelist, combined.input = combined.input, log.residual=log.residual)
#  par(kPlotParameters)
#       inds <- intersect(names(kds.1), names(kds.2))
#       kds.1 <- kds.1[inds]
#       kds.2 <- kds.2[inds]
#     xmin <- 10^(min(log10(kds.1)) - 0.1)
#     xmax <- 10^(max(log10(kds.1)) + 0.1)
#     ymin <- 10^(min(log10(kds.2)) - 0.1)
#     ymax <- 10^(max(log10(kds.2)) + 0.1)


#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))

#     xs <- xs[xs >= xmin & xs <= xmax]
#     ys <- ys[ys >= ymin & ys <= ymax]

#     # xmin <- min(xs)
#     # xmax <- max(xs)
#     # ymin <- min(ys)
#     # ymax <- max(ys)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#     yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#     if (position != FALSE) {
#       cex1 <- 1
#       pos1 <- position
#     } else {
#       cex1 <- 3.5
#       pos1 <- 1
#     }
#     plot(kds.1,
#       kds.2,
#       col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), pos1)],
#       cex = cex1,
#       pch = 19,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(ymin, ymax),
#       axes = FALSE,
#       ann = FALSE)
#     if (position == FALSE){
#       points(kds.1,  kds.2, pch = 19, cex = 2.5,
#              col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 2)])
#       points(kds.1,  kds.2, pch = 19, cex = 1.5,
#              col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 3)])
#       points(kds.1,  kds.2, pch = 19, cex = 0.5,
#              col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 4)])
#     }
#     plot_min <- max(xmin, ymin)
#     plot_max <- min(xmax, ymax)
#     segments(plot_min, plot_min, plot_max, plot_max, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=ymin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=yl,
#          labels=sapply(yl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=ys, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(mirna, font.main = 1, cex.main = 1.5, line=-2, adj=0.1)
#         cor_text <- round(
#                       cor(log(kds.1[inds]), log(kds.2[inds])),
#                       digits = 3
#                     )

#         title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)

#     title(xlab = paste0(site1," parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0(site2, " parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

# }



# PlotpairwiseFlankKdsmiRNAs <- function(mirna1, mirna2, experiment, start, stop, site, sitelist, combined.input=TRUE, log.residual=FALSE) {
#   kds.1 <- GetFlankKds(mirna1, experiment, start, stop, site, sitelist, combined.input = combined.input, log.residual=log.residual)
#   kds.2 <- GetFlankKds(mirna2, experiment, start, stop, site, sitelist, combined.input = combined.input, log.residual=log.residual)
#   par(kPlotParameters)
#       inds <- intersect(names(kds.1), names(kds.2))
#       kds.1 <- kds.1[inds]
#       kds.2 <- kds.2[inds]
#     xmin <- 10^(min(log10(kds.1)) - 0.1)
#     xmax <- 10^(max(log10(kds.1)) + 0.1)
#     ymin <- 10^(min(log10(kds.2)) - 0.1)
#     ymax <- 10^(max(log10(kds.2)) + 0.1)


#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))

#     xs <- xs[xs >= xmin & xs <= xmax]
#     ys <- ys[ys >= ymin & ys <= ymax]

#     # xmin <- min(xs)
#     # xmax <- max(xs)
#     # ymin <- min(ys)
#     # ymax <- max(ys)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#     yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#     plot(kds.1,
#       kds.2,
#       col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 1)],
#       cex = 3.5,
#       pch = 19,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(ymin, ymax),
#       axes = FALSE,
#       ann = FALSE)
#       points(kds.1,  kds.2, pch = 19, cex = 2.5, col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 2)])
#       points(kds.1,  kds.2, pch = 19, cex = 1.5, col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 3)])
#       points(kds.1,  kds.2, pch = 19, cex = 0.5, col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 4)])
#     plot_min <- max(xmin, ymin)
#     plot_max <- min(xmax, ymax)
#     segments(plot_min, plot_min, plot_max, plot_max, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=ymin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=yl,
#          labels=sapply(yl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=ys, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(site, font.main = 1, cex.main = 1.5, line=-2, adj=0.1)
#         cor_text <- round(
#                       cor(log(kds.1[inds]), log(kds.2[inds])),
#                       digits = 3
#                     )

#         title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)

#     title(xlab = paste0(mirna1," parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0(mirna2, " parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

# }



# GetAllMirnaFlankKds <- function(mirna, experiment, start, stop, sitelist,
#                           log.residual=FALSE, combined.input=TRUE, 
#                           site_list=NULL, nosite=TRUE) {
#   # Get Kds for the miRNA.
#   kds <- GetKds(mirna, experiment, start, stop, sitelist,
#                 log.residual=log.residual, combined.input=combined.input,
#                 nosite=nosite)
#   # Get the 8mer flanks, so as to populate the flank names as row names.
#   kds.8mer <- GetFlankKds(mirna, experiment, "8mer", start, stop, sitelist,
#                           log.residual=log.residual,
#                           combined.input=combined.input, nosite=nosite)
#   print(kds.8mer)
#   flanks <- names(kds.8mer)
#   # Pre-allocate the matrix with "NA".
#   flanks.all <- matrix(NA, nrow=length(kds.8mer), ncol=length(kds) - 3)
# K
#   rownames(flanks.all) <- names(kds.8mer)
#   colnames(flanks.all) <- names(kds)[1:(length(kds) - 3)]

#   # Iterate over all the kds names. Exclude the last three since this is
#   # the background parameter, the Ago concentration vector, and the
#   # loglikelihood or sum of squares for the log-transformed fit.
#   sapply(names(kds)[1:(length(kds) - 3)], function(site) {
#     kds_flanks <- try(GetFlankKds(mirna, experiment, site, start, stop, sitelist, log.residual=log.residual, combined.input=combined.input))
#     if (length(kds_flanks) > 0) {
#     flanks.all[names(kds_flanks), site] <<- as.numeric(kds_flanks)
#   }
#     })

#   flanks.new <- matrix(mapply(flanks.all, FUN=as.numeric),
#                        ncol=ncol(flanks.all), nrow=nrow(flanks.all))

#   rownames(flanks.new) <- rownames(flanks.all)
#   colnames(flanks.new) <- colnames(flanks.all)

#   return(flanks.new)
# }

# GetKineticsData <- function(mirna, experiment, n_constant, sitelist, kDilRatio, kSubSet=FALSE) {
#   sitesXcounts <- GetSitesXCounts(mirna, "equilibrium", n_constant, sitelist)
#   print(sitesXcounts)
#   sitesXcounts.kinetics <- GetSitesXCountsKinetics(mirna, experiment, n_constant,
#                                               sitelist)
#   print(sitesXcounts.kinetics)
#   sitesXcounts.p <- sitesXcounts.kinetics[[1]]
#   sitesXcounts.c <- sitesXcounts.kinetics[[2]]
#   # 2. Load the Kds from the equilibrium fits
#   params.e <- GetSiteKds(mirna, "equilibrium", n_constant, sitelist)
#   # I never use the sequences, realistically
#   seqs <- sitesXcounts[,1]
#   # Subset the equilibrium and kinetic data matrices, such that the equilibrium
#   # Matrix includes the combined input through to the zero protein, and the
#   # kinetics matrices include the pulse through to the equilibrium sample.
#   sitesXcounts <- sitesXcounts[, c(-1, -2)]
#   indeces.to.use <- which(rowSums(sitesXcounts) > 0 & sitesXcounts[, 1] > 0)
#   print(indeces.to.use)
#   kData.e <- sitesXcounts[indeces.to.use, 1 : 7]
#   kData.p <- sitesXcounts.p[indeces.to.use, 4 : (ncol(sitesXcounts.p) - 2)]
#   kData.c <- sitesXcounts.c[indeces.to.use, 4 : (ncol(sitesXcounts.c) - 2)]
#   kData.c[, 1] <- 0
#   # Conditionally subset the data if flag is not FALSE:
#   if (kSubSet != FALSE) {
#     print(kSubSet)
#     kData.e.sites <- kData.e[1:kSubSet, ]
#     kData.e.none  <- colSums(kData.e[(kSubSet + 1):nrow(kData.e), ])
#     kData.e       <- rbind(kData.e.sites, kData.e.none)
#     rownames(kData.e)[kSubSet + 1] <- "None"
#     kData.p.sites <- kData.p[1:kSubSet, ]
#     kData.p.none  <- colSums(kData.p[(kSubSet + 1):nrow(kData.p), ])
#     kData.p       <- rbind(kData.p.sites, kData.p.none)
#     rownames(kData.p)[kSubSet + 1] <- "None"
#     kData.c.sites <- kData.c[1:kSubSet, ]
#     kData.c.none  <- colSums(kData.c[(kSubSet + 1):nrow(kData.c), ])
#     kData.c       <- rbind(kData.c.sites, kData.c.none)
#     rownames(kData.c)[kSubSet + 1] <- "None"
#   }
#   if (mirna == "let-7a") {
#   colnames(kData.p)[2:5] <- c("2.1", "5.1", "2.2", "5.2")
#   colnames(kData.c) <- colnames(kData.p)

#   }
#   times <- as.integer(floor(as.numeric(colnames(kData.p))))
#   kData <- sapply(unique(times), GetAverageOfReplicates,
#                    times=times, data=rbind(kData.p, kData.c))
#   print("Kdata:")
#   print(kData)
#   print("done kData")
#   rownames(kData) <- c(sapply(c(".p", ".c"), function(i) {
#     sapply(rownames(kData.p), function(site) {
#       paste0(site, i , collapse = "")
#       })
#     }))
#   kTimes <- unique(times)/60
#   colnames(kData) <- kTimes
#   kInputPulseReads <- sitesXcounts.p[indeces.to.use, c("I_combined")]
# kInputChaseReads <- sitesXcounts.c[indeces.to.use, c("I_combined")]
# if (kSubSet != FALSE) {
#   kInputPulseReads.sites <- kInputPulseReads[1:kSubSet]
#   kInputPulseReads.none  <- sum(kInputPulseReads[(kSubSet + 1):length(kInputPulseReads)])
#   kInputPulseReads       <- c(kInputPulseReads.sites, kInputPulseReads.none)
#   kInputChaseReads.sites <- kInputChaseReads[1:kSubSet]
#   kInputChaseReads.none  <- sum(kInputChaseReads[(kSubSet + 1):length(kInputChaseReads)])
#   kInputChaseReads       <- c(kInputChaseReads.sites, kInputChaseReads.none)
# }
# # Fit the ratio of pulse to chase, in the native library,
# # to correct for the actual difference in concentration from 1:1 in the
# # experiment.
# kPulseChaseRatio <- (sum(sitesXcounts.c[indeces.to.use,c("I")]) /
#                      sum(sitesXcounts.p[indeces.to.use,c("I")]))
# kInputPulseConc <- Norm(kInputPulseReads) * kLibraryConcInRxn
# kInputChaseConc <- Norm(kInputChaseReads) * kLibraryConcInRxn * kPulseChaseRatio
# names(kInputPulseConc) <- rownames(kData.p)
# names(kInputChaseConc) <- rownames(kData.c)
# # Constants throughout the optimization routine:
# kAgoStockConc <- params.e["AGO","Mean"]
# kIndsKDs        <<- 1:kNumSites
# kIndsKoffs      <<- (kNumSites + 1):(2 * kNumSites)
# kIndsContamKD   <<- 2 * kNumSites + 1
# kIndsContamKoff <<- 2 * kNumSites + 2
# kIndsAgo        <<- 2 * kNumSites + 3
# kIndsContam     <<- 2 * kNumSites + 4
# kIndsBgs        <<- (2 * kNumSites + 5):(2 * kNumSites + 16)
# # Make the matrix to multiply through the adjoint equations:
# kFDotXMatrix <- diag(x=1,nrow=4*kNumSites+2)
# kFDotXMatrix[nrow(kFDotXMatrix) -1, 1:(2 * kNumSites)] <- 1
# kFDotXMatrix[nrow(kFDotXMatrix), (2 * kNumSites + 1):(4 * kNumSites)] <- 1
# # Define the total concentration in the INITIAL BINDING (kInputInitial), and in the
# # chase (L).
# kInputInitial <- c(kInputPulseConc, kInputChaseConc * 0)
# kInput <- ((kInputInitial + c(kInputPulseConc * 0, kInputChaseConc) * kDilRatio)
#            / (kDilRatio + 1))
# kInputMatrix <- matrix(c(kInputInitial / (1 + kDilRatio), 
#                          rep(kInput, kNumBgs - 1)),
#                        nrow=length(kInput),
#                        ncol=kNumBgs,
#                        byrow=FALSE)
# rownames(kInputMatrix) <- rownames(kData)
# colnames(kInputMatrix) <- colnames(kData)
# kInputTotals <- colSums(kInputMatrix)
#   return(list(kData,kInputMatrix))
# }

# GetBgs <- function(pars, kNumSites) {
#   kIndsBgs        <- (2 * kNumSites + 5):(2 * kNumSites + 16)
#   return(10^pars[kIndsBgs])
# }

# MakeODEPars <- function(pars, kNumSites, l. = kInput) {
#   pars <- 10^pars
#   kIndsKDs        <- 1:kNumSites
#   kIndsKoffs      <- (kNumSites + 1):(2 * kNumSites)
#   kIndsContamKD   <- 2 * kNumSites + 1
#   kIndsContamKoff <- 2 * kNumSites + 2
#   # The site type on and off rates, the background terms for each column,
#   # the and the total concentration of Ago (a) and the contaminant (b):
#   print(kIndsKDs)
#   print(kIndsKoffs)
#   kds.a.   <- pars[kIndsKDs]
#   koffs.a. <- pars[kIndsKoffs]
#   kd.b     <- pars[kIndsContamKD]
#   koff.b   <- pars[kIndsContamKoff]
#   # Assignment of the Kds (being the ratio of the on and off rates)
#   kons.a. <- koffs.a. / kds.a.
#   kon.b   <- koff.b / kd.b
#   # Removal of Kd parameter as it is no longer useful:
#   # Calculate the vector of diluted, bound RNA:
#   # Assign the parameter vector and initial conditions vector for the ODE solver:
#   parms.k <- c(koffs.a., kons.a., l., koff.b, kon.b)
#   return(parms.k)
# }

# MakeX0 <-function(pars, kNumSites, l0. = kInputInitial) {
#   pars <- 10^pars
#   kIndsKDs        <- 1:kNumSites
#   kIndsKoffs      <- (kNumSites + 1):(2 * kNumSites)
#   kIndsContamKD   <- 2 * kNumSites + 1
#   kIndsContamKoff <- 2 * kNumSites + 2
#   kIndsAgo        <- 2 * kNumSites + 3
#   kIndsContam     <- 2 * kNumSites + 4

#   # The site type on and off rates, the background terms for each column,
#   # the and the total concentration of Ago (a) and the contaminant (b):
#   kds.a.   <- pars[kIndsKDs]
#   koffs.a. <- pars[kIndsKoffs]
#   kd.b     <- pars[kIndsContamKD]
#   koff.b   <- pars[kIndsContamKoff]
#   A0       <- 0.4 * pars[kIndsAgo]
#   B0       <- 0.4 * pars[kIndsContam]
#   # Assignment of the Kds (being the ratio of the on and off rates)
#   kons.a. <- koffs.a. / kds.a.
#   kon.b   <- koff.b / kd.b
#   # print(kons.a.)
#   # print("kons.a.")
#   # print(kon.b)
#   # print("kon.b")
#   # print(kds.a.)
#   # print("kds.a.")
#   # print(koffs.a.)
#   # print("koffs.a.")
#   # print(kd.b)
#   # print("kd.b")
#   # print(A0)
#   # print("A0")
#   # print(B0)
#   # print("B0")
#   # print(l0.)
#   # print("l0.")
#   # Get the free amount of Ago and contaminant:
#   a0.b0. <- GetFreeAgoAndContam(kds.a., kd.b, l0., A0, B0)
#   # print("a0.b0.")
#   # print(a0.b0.)
#   a0 <- a0.b0.[1]
#   b0 <- a0.b0.[2]
#   a0 <<- a0
#   b0 <<- b0
#   l0. <<- l0.
#   # print("a0")
#   # print(a0)
#   # print("b0")
#   # print(b0)

#   # Calculate the initial Ago and contaminant occupancies:
#   a.occs0.b.occs0 <- GetOccupanciesContaminant(a0, b0, kds.a., kd.b)
#   xap0.xac0. <- l0. * a.occs0.b.occs0[[1]]
#   xbp0.xbc0. <- l0. * a.occs0.b.occs0[[2]]
#   # print("xap0.xac0.")
#   # print(xap0.xac0.)
#   # print("xbp0.xbc0.")
#   # print(xbp0.xbc0.)
#   # Removal of Kd parameter as it is no longer useful:
#   # Calculate the vector of diluted, bound RNA:
#   # Assign the parameter vector and initial conditions vector for the ODE solver:
#   x0. <- c(xap0.xac0., xbp0.xbc0., a0.b0.) / (kDilRatio + 1)
#   return(x0.)
# }

# MakeXTimeCourse <- function(x0, parms.k, kNumSites, kCScriptDir, kCScriptName, kTimesAll, kTimes, verbose.=FALSE) {
#   # Load and run the ODE:
#   dyn.load(kCScriptDir)
#   # print(x0)
#   # print(parms.k)
#   x.t <- t(ode(y          = x0,
#                times      = kTimesAll,
#                func       = "ode_deriv",
#                parms      = parms.k,
#                method     = "lsodes",
#                dllname    = kCScriptName,
#                sparsetype = "sparseint",
#                initfunc   = "ode_p_init",
#                nout       = 1,
#                verbose    = verbose.)[,2:(4*kNumSites+3)])
#   dyn.unload(kCScriptDir)
#   x <- x.t[, which(kTimesAll %in% kTimes)]
#   print(x)
#   print(kTimes)
#   colnames(x) <- kTimes
#   # Calculate total bound pulse and chase sites:
#   return(x)
# }


# MakeModelPrediction <- function(x, bgs, data = kData, l = kInputMatrix) {
#   # Combine the ago-bound and contaminant-bound site types:
#   # Written out the numerator for maximum accuracy:
#   # Form of equation is :        x(L - X) + B(l - x)
#   #                          D * –––––––––––––––––––
#   #                                (L - X)(X + B)
#   kNumSites = nrow(kData)/2
#   # print(dim(x))
#   # print(dim(l))
#   # print(length(bgs))
#   # print(dim(data))
#   x.ago <- x[1:(2 * kNumSites), ] 
#   x.con <- x[(2 * kNumSites + 1):(4 * kNumSites), ]
#   x     <- x.ago + x.con
#   X     <- colSums(x)
#   L     <- colSums(l)
#   D     <- colSums(kData)
#   B     <- bgs
#   B <<- B
#   L <<- L
#   X <<- X
#   # Transpose D-multiply x and l for row-multiplication:
#   Dx.     <- D * t(x)
#   Dl.     <- D * t(l)
#   Dl. <<- Dl.
#   Dx. <<- Dx.
#   # print(Dx.)
#   # print(Dl.)
#   # print(B)
#   pred.  <- (Dx.*L - Dx.*X + Dl.*B - Dx.*B) /
#                 (L*B + L*X -  B*X - X^2)
#   pred   <- t(pred.)
#   return(pred)
# }


# GetKineticsModelGeneral <- function(mirna, experiment, n_constant, sitelist, costfunction,dil=FALSE,pars=FALSE) {
#   if (pars == FALSE) {
#     pars <- GetSiteKoffs(mirna, experiment, n_constant, sitelist, costfunction, dil=dil)
#   }
#   if (dil == TRUE) {
#     kDilRatio <- 10^pars[length(pars)-12]
#     print(kDilRatio)
#   } else {
#     kDilRatio <- 11
#   }
#   kData_kInput <- GetKineticsData(mirna, experiment, n_constant, sitelist, kDilRatio=11)
#   kData <- kData_kInput[[1]]
#   kNumSites <- nrow(kData)/2
#   print(kNumSites)
#   # system(paste0("python SolveForOffRates/",
#   #             "MakeCScriptSingleExponential.py ",
#   #             kNumSites))
#   # # Performa a system command to compile both of these scripts.
#   # system(paste0("R CMD SHLIB SolveForOffRates/SingleODEContam_",
#   #               kNumSites, ".c"))

#   kInputMatrix <- kData_kInput[[2]]

# # Assign names to their name and directory to be used within the optimization
# # routine.
#   kCScriptName <- paste0("SingleODEContam_", kNumSites)
#   kCScriptDir <- paste0("SolveForOffRates/", kCScriptName, ".so")
#   try(dyn.unload(kCScriptDir))
#   print("hi")
#   print(names(pars))
#   temp_names <- names(pars)
#   pars <- as.numeric(pars)

#   names(pars) <- temp_names
#   if (dil != FALSE) {
#     kDilRatio <- 10^pars[2 * kNumSites + 5]
#     pars <- pars[-(2*kNumSites + 5)]
#   } else {
#     kDilRatio <- 11
#   }
#   x0      <- MakeX0(pars, kNumSites, l0. = (kDilRatio+1) * kInputMatrix[,1])
#   x0 <<- x0
#   print("done x0")
#   # print(x0)
#   parms.k <- MakeODEPars(pars, kNumSites, l. = kInputMatrix[, 2])
#   print("done parms.k")
#   parms.k <<- parms.k
#   # print(parms.k)
#   x       <- MakeXTimeCourse(x0, parms.k, kNumSites, kCScriptDir, kCScriptName,kTimesAll, kTimes)
#   print("done x")
#   x <<- x
#   kData <<- kData
#   bgs     <- GetBgs(pars, kNumSites)
#   bgs <<- bgs
#   model     <- MakeModelPrediction(x, bgs, data = kData, l = kInputMatrix)
#   print("done model")
#   return(model)
# }

# PlotKineticsFit <- function(mirna, num_sites, sitelist, costfunc, dil=FALSE, col=FALSE) {
#   kinetics.data <- GetKineticsData(mirna, "kinetics", num_sites, sitelist,kDilRatio=11)
#   kData <- kinetics.data[[1]]
#   kInputMatrix <- kinetics.data[[2]]
#   data.e <- GetSitesXCounts(mirna, "equilibrium", 5, "paper")
#   colors <- kSiteColors[rownames(data.e), ]
#   model <- GetKineticsModelGeneral(mirna, "kinetics", 5, "paper", costfunc, dil=dil)
#   model <<- model
#   par(kPlotParameters)
#   model.norm <- t(t(model + 1) / colSums(model + 1))
#   data.norm <- t(t(kData + 1) / colSums(kData + 1))
#   if (col != FALSE) {
#     model.norm <- model.norm[, col, drop=FALSE]
#     data.norm <- data.norm[, col, drop=FALSE]
#   }
#   plot(c(model.norm), c(data.norm),
#        xlim = c(1e-8, 0.8),
#        ylim = c(1e-8, 0.8),
#        log  = 'xy',
#        col  = colors,
#        pch  = rep(c(20,1),
#        each = kNumSites),
#        cex = 1) 
#     abline(0, 1, lty = 3)
# }

# PlotKineticsFitInidividual <- function(mirna, n_constant, sitelist, costfunc, log. = "", dil=FALSE) {
#   site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", sitelist,".txt"),
#     stringsAsFactors=FALSE)[,1], "None")
#   print(site_list)
#   kinetics.data <- GetKineticsData(mirna, "kinetics", n_constant, sitelist, kDilRatio=11)
#   kData <- kinetics.data[[1]]
#   kInputMatrix <- kinetics.data[[2]]
#   data.e <- GetSitesXCounts(mirna, "equilibrium", 5, sitelist)
#   model <- GetKineticsModelGeneral(mirna, "kinetics", 5, sitelist, costfunc, dil = dil)
#   model <<- model
#   koffs <- GetSiteKoffs(mirna, "kinetics", n_constant, sitelist, costfunc, dil=dil)
#   print(koffs)
#   koffs <- 10^koffs[(nrow(kData)/2+1):(nrow(kData))]
#   names(koffs) <- rownames(data.e)
#   print(koffs)

#   times <- as.numeric(colnames(model))
#   print(times)
#   if (sitelist == "canonical") {
#     dev.new(xpos = 20, ypos = 20, height = 7, width = 10)
#     par(mfrow = c(2, 4))
#   } else if (mirna == "miR-124") {

#     dev.new(xpos = 20, ypos = 20, height = 8, width = 14)
#     par(mfrow = c(4, 6))  
#   } else {
#     dev.new(xpos = 20, ypos = 20, height = 8, width = 12)
#     par(mfrow = c(4, 5))  
#   }
#   par(kPlotParameters)

#   kInput.Pulse <- Norm(kInputMatrix[1:(nrow(kData)/2),1])
#   kInput.Chase <- Norm(kInputMatrix[(nrow(kData)/2 + 1):nrow(kData),2])

#   kInput.Pulse <<- kInput.Pulse
#   kInput.Chase <<- kInput.Chase
#   kInput.Pulse <- 1
#   kInput.Chase <- 1
#   for (site in site_list) {
#     print(site)
#     color <- kSiteColors[site, ]
#     times[1] <- 0.5/60
#     site.inds <- grep(paste0(site,"."), rownames(model),fixed = TRUE)
#     print(rownames(model))
#     print(site.inds)
#     model.p <- model[site.inds[1], ]/colSums(model)
#     model.c <- model[site.inds[2], ]/colSums(model)
#     data.p <- kData[site.inds[1], ]/colSums(kData)
#     data.c <- kData[site.inds[2], ]/colSums(kData)
#     # model.c[1] <- 1
#     # data.c[1] <- 1

#     print(model.p)
#     print(model.c)
#     print(data.p)
#     print(data.c)

#     ylim. = c(min(model.p, model.c[-1], data.p, data.c[-1]),
#              max(model.p, model.c, data.p, data.c))
#     # ylim. <- c(0.1, max(model.p, model.c, data.p, data.c)*2)
#     print(ylim.)
#     plot(times, data.p,
#            xlim = c(0.3, 20000)/60,
#            ylim = ylim.,
#            log  = log.,
#            col  = color,
#            pch  = 20) 
#     points(times, data.c, pch = 1,col=color)
#     lines(times, model.p, col=color)
#     lines(times, model.c, lty = 2, col=color)
#     mtext(site, 3, line = 0.5, adj=0.4)
#     print(grep(paste0(site),names(koffs),fixed=TRUE))
#     dwell = 1/(koffs[site])
#     segments(dwell,ylim.[1],dwell,ylim.[2], lty=2,col=color)
#   }
#     plot(1, type = "n", ann=FALSE, axes=FALSE)
#   title(main = costfunc, cex = 2)

# }
# # graphics.off()
# # PlotKoffsVsKds("miR-1", 5, "paper", dil=TRUE)
# # dev.copy2pdf(file = "171208_miR-1_koffs_logres.pdf")
# # graphics.off()
# # PlotKoffsVsKds("let-7a", 5, "paper", dil=TRUE)
# # dev.copy2pdf(file = "171208_let-7a_koffs_logres.pdf")
# # graphics.off()
# # PlotKoffsVsKds("miR-124", 5, "paper", dil=TRUE)
# # dev.copy2pdf(file = "171208_miR-124_koffs_logres.pdf")
# # graphics.off()
# # PlotKoffsVsKds("lsy-6", 5, "paper", dil=TRUE)
# # dev.copy2pdf(file = "171208_lsy-6_koffs_logres.pdf")
# # graphics.off()
# # PlotKineticsFitInidividual("miR-1", 5, "paper", "logres", log.="xy", dil=TRUE)
# # dev.copy2pdf(file = "171208_miR-1_fits_logres.pdf")
# # graphics.off()
# # PlotKineticsFitInidividual("let-7a", 5, "paper", "logres", log.="xy", dil=TRUE)
# # dev.copy2pdf(file = "171208_let-7a_fits_logres.pdf")
# # graphics.off()
# # PlotKineticsFitInidividual("miR-124", 5, "paper", "logres", log.="xy", dil=TRUE)
# # dev.copy2pdf(file = "171208_miR-124_fits_logres.pdf")
# # graphics.off()
# # PlotKineticsFitInidividual("lsy-6", 5, "paper", "logres", log.="xy", dil=TRUE)
# # dev.copy2pdf(file = "171208_lsy-6_fits_logres.pdf")
# # graphics.off()

# # PlotKineticsFitInidividual("miR-1", 5, "paper", "multinom", log.="xy", dil=TRUE)
# # dev.copy2pdf(file = "171208_miR-1_fits_multinom.pdf")
# # graphics.off()
# # PlotKineticsFitInidividual("let-7a", 5, "paper", "multinom", log.="xy", dil=TRUE)
# # dev.copy2pdf(file = "171208_let-7a_fits_multinom.pdf")
# # graphics.off()
# # PlotKineticsFitInidividual("miR-124", 5, "paper", "multinom", log.="xy", dil=TRUE)
# # dev.copy2pdf(file = "171208_miR-124_fits_multinom.pdf")
# # graphics.off()
# # PlotKineticsFitInidividual("lsy-6", 5, "paper", "multinom", log.="xy", dil=TRUE)
# # dev.copy2pdf(file = "171208_lsy-6_fits_multinom.pdf")
# # graphics.off()



# PlotKineticsFitInidividualPars <- function(mirna, n_constant, sitelist, costfunc, log. = "", dil=FALSE,pars=FALSE) {
#   site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", sitelist,".txt"),
#     stringsAsFactors=FALSE)[,1], "None")
#   print(site_list)
#   kinetics.data <- GetKineticsData(mirna, "kinetics", n_constant, sitelist,kDilRatio=11)
#   kData <- kinetics.data[[1]]
#   kInputMatrix <- kinetics.data[[2]]
#   data.e <- GetSitesXCounts(mirna, "equilibrium", 5, sitelist)
#   model <- GetKineticsModelGeneral(mirna, "kinetics", 5, sitelist, costfunc, dil = dil, pars=pars)
#   print(model)
#   print("done model")
#   koffs <- GetSiteKoffs(mirna, "kinetics", n_constant, sitelist, costfunc, dil=dil)
#   print(koffs)
#   koffs <- 10^koffs[(nrow(kData)/2+1):(nrow(kData))]
#   names(koffs) <- rownames(data.e)
#   print(koffs)

#   times <- colnames(model)
#   if (sitelist == "canonical") {
#     dev.new(xpos = 20, ypos = 20, height = 7, width = 10)
#     par(mfrow = c(2, 4))
#   } else {
#     dev.new(xpos = 20, ypos = 20, height = 10, width = 17)
#     par(mfrow = c(4, 6))  
#   }
#   par(kPlotParameters)

#   for (site in site_list) {
#     print(site)
#     color <- kSiteColors[site, ]
#     times[1] <- 0.5/60
#     site.inds <- grep(paste0(site,"."), rownames(model),fixed = TRUE)
#     print(rownames(model))
#     print(site.inds)
#     model.p <- model[site.inds[1], ]/colSums(model)
#     model.c <- model[site.inds[2], ]/colSums(model)
#     data.p <- kData[site.inds[1], ]/colSums(kData)
#     data.c <- kData[site.inds[2], ]/colSums(kData)
#     # print(model.p)
#     # print(model.c)
#     # print(data.p)
#     # print(data.c)
#     ylim. = c(min(model.p, model.c[-1], data.p, data.c[-1]),
#              max(model.p, model.c, data.p, data.c))
#     print(ylim.)

#     plot(times, data.p,
#            xlim = c(0.3, 20000)/60,
#            ylim = ylim.,
#            log  = log.,
#            col  = color,
#            pch  = 20) 
#     points(times, data.c, pch = 1,col=color)
#     lines(times, model.p, col=color)
#     lines(times, model.c, lty = 2, col=color)
#     mtext(site, 3, line = 2, adj=0.4)
#     print(grep(paste0(site),names(koffs),fixed=TRUE))
#     dwell = 1/(koffs[site])
#     segments(dwell,ylim.[1],dwell,ylim.[2], lty=2,col=color)
#   }
# }




# MakeFigure1 <- function() {
#   graphics.off()
#   PlotSiteScatterWithInput("miR-1", "equilibrium", 5, "canonical", 7, xpos = 20, ypos = 20)
#   dev.copy2pdf(file = "2017_Paper/Figure1B_raw_v6.eps")
#   PlotSiteEnrichments("miR-1", "equilibrium", 5, "canonical", "canonical", xpos = 620, ypos = 20)
#   dev.copy2pdf(file = "2017_Paper/Figure1C_raw_v6.eps")
#   PlotSiteOccupancy("miR-1", "equilibrium", 5, "canonical", "canonical", xpos = 20, ypos = 420)
#   dev.copy2pdf(file = "2017_Paper/Figure1D_raw_v6.eps")
#   PlotSiteEnrichments("miR-1", "equilibrium", 5, "paper", "paper", xpos = 620, ypos = 420)
#   dev.copy2pdf(file = "2017_Paper/Figure1E_raw_v6.eps")
#   PlotSiteKds("miR-1", "equilibrium", 5, "paper", "paper", xpos = 20, ypos = 820)
#   dev.copy2pdf(file = "2017_Paper/Figure1F_raw_v6.eps")
# }

# MakeFigure2 <- function() {
#   graphics.off()
#   PlotSiteKds("let-7a", "equilibrium", 5, "paper", "paper", xpos = 20, ypos = 20)
#   dev.copy2pdf(file = "2017_Paper/Figure2A_raw_v6.eps")
#   PlotSiteKds("miR-155", "equilibrium", 5, "paper", "paper", xpos = 820, ypos = 20)
#   dev.copy2pdf(file = "2017_Paper/Figure2B_raw_v6.eps")
#   PlotSiteKds("miR-124", "equilibrium", 5, "paper", "paper", xpos = 20, ypos = 620)
#   dev.copy2pdf(file = "2017_Paper/Figure2C_raw_v6.eps")
#   PlotSiteKds("lsy-6", "equilibrium", 5, "paper", "paper", xpos = 820, ypos = 620)
#   dev.copy2pdf(file = "2017_Paper/Figure2D_raw_v6.eps")
# }

# MakeSupplementalFigure2 <- function() {
#   graphics.off()
#   PlotBaekKds("miR-1", "equilibrium", 5, xpos = 20, ypos = 20)
#   dev.copy2pdf(file = "2017_Paper/FigureS2A_raw_v6.eps")
#   PlotBaekKds("let-7a", "equilibrium", 5, xpos = 820, ypos = 20)
#   dev.copy2pdf(file = "2017_Paper/FigureS2B_raw_v6.eps")
#   PlotBaekKds("miR-155", "equilibrium", 5, xpos = 20, ypos = 220)
#   dev.copy2pdf(file = "2017_Paper/FigureS2C_raw_v6.eps")
#   PlotBaekKds("miR-124", "equilibrium", 5, xpos = 820, ypos = 220)
#   dev.copy2pdf(file = "2017_Paper/FigureS2D_raw_v6.eps")
#   PlotBaekKds("lsy-6", "equilibrium", 5, xpos = 20, ypos = 420)
#   dev.copy2pdf(file = "2017_Paper/FigureS2E_raw_v6.eps")
#   PlotPositionalKds("equilibrium", 5, "centered11", xpos = 820, 420)
#   dev.copy2pdf(file = "2017_Paper/FigureS2F_raw_v6.eps")
# }

# MakeFigure3 <- function() {
#   graphics.off()
#   PlotSiteFlankEnrichments("miR-1", "equilibrium", 5, "paper", "paper","8mer", xpos = 20, ypos = 20)
#   dev.copy2pdf(file = "2017_Paper/Figure3A_raw_v7.eps")
#   PlotSiteFlankKds("miR-1", "equilibrium", 5, "paper", "paper", xpos = 820, ypos = 20)
#   dev.copy2pdf(file = "2017_Paper/Figure3B_raw_v7.eps")
#   PlotCanonicalSiteFlankingKdHeatmap("equilibrium", 5, "paper", xpos = 20, ypos = 420)
#   dev.copy2pdf(file = "2017_Paper/Figure3C_raw_v7.eps")
#   PlotStructureVsFlankingKds("miR-1", "equilibrium", "I_combined", 5, "paper", "8mer", 1, 15, xpos=20, ypos=20, absolute=TRUE,noconstant=FALSE)
#   dev.copy2pdf(file = "2017_Paper/Figure3D_raw_v7.eps")
#   PlotPlSampleFlanks("miR-1", "equilibrium", 0.4, 5, "paper", "8mer", 1, 15, xpos=520, ypos=20)
#   dev.copy2pdf(file = "2017_Paper/Figure3E_raw_v7.eps")
# }

# MakeFigure4 <- function(mirna,site) {
#   graphics.off()
#   # ypos = 20
#   # PlotSiteKdsVsRepression("miR-1","equilibrium", 5, "paper", noncanon=FALSE,cutoff=FALSE,          xpos =  20, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4A1_raw_v8.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("let-7a","equilibrium", 5, "paper", noncanon=FALSE,cutoff=FALSE,          xpos =  420, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4A2_raw_v8.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-155","equilibrium", 5, "paper", noncanon=FALSE,cutoff=FALSE,          xpos =  820, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4A3_raw_v8.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-124","equilibrium", 5, "paper", noncanon=FALSE,cutoff=FALSE,          xpos =  1220, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4A4_raw_v8.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("lsy-6","equilibrium", 5, "paper", noncanon=FALSE,cutoff=FALSE,          xpos =  1620, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4A5_raw_v8.pdf", useDingbats=FALSE)
  
#   # cutoff. <- 0
#   # PlotSiteKdsVsRepression("miR-1","equilibrium", 5, "paper", noncanon=TRUE,cutoff=FALSE,          xpos =  20, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b1_raw_v8_alt1.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("let-7a","equilibrium", 5, "paper", noncanon=TRUE,cutoff=FALSE,          xpos =  420, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b2_raw_v8_alt1.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-155","equilibrium", 5, "paper", noncanon=TRUE,cutoff=FALSE,          xpos =  820, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b3_raw_v8_alt1.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-124","equilibrium", 5, "paper", noncanon=TRUE,cutoff=FALSE,          xpos =  1220, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b4_raw_v8_alt1.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("lsy-6","equilibrium", 5, "paper", noncanon=TRUE,cutoff=FALSE,          xpos =  1620, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b5_raw_v8_alt1.pdf", useDingbats=FALSE)

#   # cutoff. <- -1
#   # PlotSiteKdsVsRepression("miR-1","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  20, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b1_raw_v8_alt2.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("let-7a","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  420, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b2_raw_v8_alt2.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-155","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  820, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b3_raw_v8_alt2.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-124","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  1220, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b4_raw_v8_alt2.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("lsy-6","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  1620, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b5_raw_v8_alt2.pdf", useDingbats=FALSE)
#   # ypos = 820
#   # cutoff. <- -2
#   # PlotSiteKdsVsRepression("miR-1","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  20, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b1_raw_v8_alt3.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("let-7a","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  420, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b2_raw_v8_alt3.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-155","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  820, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b3_raw_v8_alt3.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-124","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  1220, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b4_raw_v8_alt3.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("lsy-6","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  1620, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b5_raw_v8_alt3.pdf", useDingbats=FALSE)

#   # cutoff. <- -3
#   # PlotSiteKdsVsRepression("miR-1","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  20, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b1_raw_v8_alt4.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("let-7a","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  420, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b2_raw_v8_alt4.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-155","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  820, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b3_raw_v8_alt4.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-124","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  1220, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b4_raw_v8_alt4.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("lsy-6","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  1620, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b5_raw_v8_alt4.pdf", useDingbats=FALSE)

#   # cutoff. <- -4
#   # PlotSiteKdsVsRepression("miR-1","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  20, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b1_raw_v8_alt5.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("let-7a","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  420, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b2_raw_v8_alt5.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-155","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  820, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b3_raw_v8_alt5.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("miR-124","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  1220, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b4_raw_v8_alt5.pdf", useDingbats=FALSE)
#   # PlotSiteKdsVsRepression("lsy-6","equilibrium", 5, "paper", noncanon=TRUE,cutoff=cutoff.,          xpos =  1620, ypos = ypos)
#   # dev.copy2pdf(file = "2017_Paper/Figure4b5_raw_v8_alt5.pdf", useDingbats=FALSE)


#   # PlotAllSiteKdsVsRepression("equilibrium", 5, "paper", noncanon=FALSE,cutoff=FALSE,                xpos = 2620, ypos = 520)
#   # dev.copy2pdf(file = "2017_Paper/Figure4B1_raw_v8.pdf", useDingbats=FALSE)
#   # PlotAllSiteKdsVsRepression("equilibrium", 5, "paper", noncanon=TRUE,cutoff=FALSE,                xpos = 2620, ypos = 520)
#   # dev.copy2pdf(file = "2017_Paper/Figure4B2_raw_v8.pdf", useDingbats=FALSE)
#   # PlotAllSiteKdsVsRepression("equilibrium", 5, "paper", noncanon=TRUE,cutoff=-1,                xpos = 2620, ypos = 520)
#   # dev.copy2pdf(file = "2017_Paper/Figure4B3_raw_v8.pdf", useDingbats=FALSE)
#   # PlotAllSiteKdsVsRepression("equilibrium", 5, "paper", noncanon=TRUE,cutoff=-2,                xpos = 2620, ypos = 520)
#   # dev.copy2pdf(file = "2017_Paper/Figure4B4_raw_v8.pdf", useDingbats=FALSE)
#   # PlotAllSiteKdsVsRepression("equilibrium", 5, "paper", noncanon=TRUE,cutoff=-3,                xpos = 2620, ypos = 520)
#   # dev.copy2pdf(file = "2017_Paper/Figure4B5_raw_v8.pdf", useDingbats=FALSE)
#   # PlotAllSiteKdsVsRepression("equilibrium", 5, "paper", noncanon=TRUE,cutoff=-4,                xpos = 2620, ypos = 520)
#   # dev.copy2pdf(file = "2017_Paper/Figure4B6_raw_v8.pdf", useDingbats=FALSE)


#   # num_tick <- 1
#   # for (mirna in mirnas.all) {
#   #   for (site in canonical.sites) {
#   #     PlotSingleFlankKdsVsRepression(mirna, "equilibrium", 5, "paper", site)
#   #     print(paste0("2017_Paper/Figure4C", num_tick, "_raw_v8.pdf"))
#   #     dev.copy2pdf(file = paste0("2017_Paper/Figure4C", num_tick, "_raw_v8.pdf"), useDingbats=FALSE)
#   #     num_tick <- num_tick + 1
  
#   #   }
#   # }

#   num_tick <- 1
#   for (mirna in mirnas.all) {
#     PlotAllCanonicalKdsVsRepression(mirna, "equilibrium", 5, "paper")
#     dev.copy2pdf(file = paste0("2017_Paper/Figure4D", num_tick, "_raw_v8.pdf"), useDingbats=FALSE)    
#     num_tick <- num_tick + 1
#   }
#   # PlotLinearModelKdVsRepressionCanonical("equilibrium", 5, "paper")
# }

# #   k.R <- 1.987e-3 # in kcal K-1 mol-1
# #   k.T <- 310
# # dG.df <- read.table("canonical.sites_mfe.txt")
# # dG.df.new <- read.table("canonical.sites_mfe_new.txt")

# # canon_2.7 <- c(2, 4, 6)

# # plotkdsvsEnergy <- function(xpos = 20, ypos = 20){
# #   # graphics.off()
# #   dev.new(height = 10, width = 4, xpos = xpos , ypos = ypos)
# #   par(kPlotParameters)
# #   par(mfrow=c(3, 1))
# #   x <- 10^seq(-12, 1, length = 10)
# #   GetConstant <- function(C, kds, del_g) {
# #     del_g_model <- k.R*k.T*log(kds[canon_2.7]) + C
# #     return(sum((del_g_model - del_g[canon_2.7])^2))
# #   }
# #     total_x <- c()
# #     total_y <- c()

# #   for (mirna in mirnas.all) {
# #     print(mirna)
# #     print(canonical.sites)
# #     kds <- GetSiteKds(mirna, "equilibrium", 5, "paper")[canonical.sites,]$Mean*10^-9
# #     mirna_new <- paste(unlist(strsplit(mirna,"-")),collapse=".")
# #     print(mirna_new)
# #     del_g <- dG.df[canonical.sites,mirna_new]
# #     total_x <- c(total_x, log(kds)[canon_2.7])
# #     total_y <- c(total_y, del_g[canon_2.7])
# #     if (mirna == "miR-1") {
# #       plot(kds[canon_2.7]*1e9, del_g[canon_2.7],col=kMirnaColors[mirna],log='x', xlim = c(1e-3, 2), ylim = c(-10, -2), pch = c(1, 2, 3))    
# #     } else {
# #       points(kds[canon_2.7]*1e9, del_g[canon_2.7],col=kMirnaColors[mirna], pch = c(1, 2, 3))    

# #     }
# #   }
# #   text(0.002, -4.5, round(cor(total_x, total_y)^2, 2))
# #   total_x <- c()
# #   total_y <- c()

# #   title(xlab = "Kd ",line=1)
# #   title(ylab = "∆G, predicted", line=1)
# #   legend("bottomright", legend = mirnas.all, col=kMirnaColors[mirnas.all], bty="n",pch=19)
# #   legend("topleft", legend = c("7mer-m8", "6mer", "6mer-m8"), bty="n",pch=c(1, 2, 3))


# #   for (mirna in mirnas.all) {
# #     print(canonical.sites)
# #     kds <- GetSiteKds(mirna, "equilibrium", 5, "paper")[canonical.sites,]$Mean*10^-9
# #     kds_bg <- GetSiteKds(mirna, "equilibrium", 5, "paper")["None",]$Mean*10^-9
# #     mirna_new <- paste(unlist(strsplit(mirna,"-")),collapse=".")
# #     del_g <- dG.df[canonical.sites,mirna_new]
# #     kds_m8 <- kds[c(1, 2, 3)]
# #     kds_nom8 <- kds[c(4, 4, 4)]
# #     del_g_rel <- del_g[c(1, 2, 3)]-del_g[c(4, 4, 4)]
# #     print(del_g_rel)
# #     kds_rel <- kds_nom8/kds_m8
# #     print(kds_rel)
# #     total_x <- c(total_x, log(kds_rel))
# #     total_y <- c(total_y, del_g_rel)
# #     # b <- optimize(GetConstant,c(-10,20),kds=kds_rel, del_g=del_g)$minimum
# #     # print(b)
# #     # y <- k.R*k.T*log(x) + b
# #     if (mirna == "miR-1") {
# #       plot(kds_rel, del_g_rel, col=kMirnaColors[mirna],log='x',  xlim = c(1, 2000), ylim=c(0,-5),pch = c(0, 1, 2))    
# #     } else {
# #       points(kds_rel, del_g_rel, col=kMirnaColors[mirna], pch = c(0, 1, 2))    

# #     }

# #     # title(main=mirna)
# #     # lines(x,y)
# #     # plot(kds_rel, del_g,col=kSiteColors[canonical.sites,],ylim = c(-15, 5), log='x')
# #   # lines(x,y)
# # }
# #   text(3, -3.5, round(cor(total_x, total_y)^2, 3))
# #   total_x <- c()
# #   total_y <- c()

# #     title(xlab = "Kd rel", line=1)
# #     title(ylab = "∆∆G, predicted",line=1)
# #   legend("bottomright", legend = mirnas.all, col=kMirnaColors[mirnas.all], bty="n",pch=19)
# #   legend("topleft", legend = c("8mer vs 6mer", "7mer-m8 vs 6mer", "7mer-A1 vs 6mer"), bty="n",pch=c(0, 1, 2))


# #     x <- 10^(seq(0,3, length = 10))
# #     # print(mean(kds_noA/kds_A))
# #     y <- -k.R * k.T *log(x)
# #       lines(x, y, lty = 2, col = "gray")

# #   for (mirna in mirnas.all) {
# #     print(canonical.sites)
# #     kds <- GetSiteKds(mirna, "equilibrium", 5, "paper")[canonical.sites,]$Mean*10^-9
# #     kds_bg <- GetSiteKds(mirna, "equilibrium", 5, "paper")["None",]$Mean*10^-9
# #     mirna_new <- paste(unlist(strsplit(mirna,"-")),collapse=".")
# #     del_g <- dG.df[canonical.sites,mirna_new]
# #     kds_rel8 <- kds[c(4)]/kds[c(1)]
# #     kds_relArelm8 <- (kds[c(4)]^2)/(kds[c(2)]*kds[c(3)])
# #     print(del_g_rel)
# #     kds_rel <- kds_nom8/kds_m8
# #     print(kds_rel)
# #     total_x <- c(total_x, log(kds_rel8))
# #     total_y <- c(total_y, log(kds_relArelm8))
# #     # b <- optimize(GetConstant,c(-10,20),kds=kds_rel, del_g=del_g)$minimum
# #     # print(b)
# #     # y <- k.R*k.T*log(x) + b
# #     if (mirna == "miR-1") {
# #       plot(kds_rel8, kds_relArelm8, col=kMirnaColors[mirna],log='xy',  xlim = c(10, 200), ylim=c(10,100),pch = c(1, 2))    
# #     } else {
# #       points(kds_rel8, kds_relArelm8, col=kMirnaColors[mirna], pch = c(1, 2))    

# #     }

# #     # title(main=mirna)
# #     # plot(kds_rel, del_g,col=kSiteColors[canonical.sites,],ylim = c(-15, 5), log='x')
# #   # lines(x,y)
# # }
# #     title(xlab = "Predicted 8mer/6mer Kd rel 7mer-A1 & 7mer-m8", line=1)
# #     title(ylab = "Actual 8mer/6mer Kd rel from 7mer-A1 & 7mer-m8", line=1.5)
# #   text(15, 90, round(cor(total_x, total_y)^2, 3))

# #     x <- 10^(seq(0,3, length = 10))
# #     # print(mean(kds_noA/kds_A))
# #     y <- -k.R * k.T *log(x)
# #       lines(x, x, lty = 2, col = "gray")




# #   legend("bottomright", legend = mirnas.all, col=kMirnaColors[mirnas.all], bty="n",pch=19)

# # #   for (mirna in mirnas.all) {
# # #     print(canonical.sites)
# # #     kds <- GetSiteKds(mirna, "equilibrium", 5, "paper")[canonical.sites,]$Mean*10^-9
# # #     kds_bg <- GetSiteKds(mirna, "equilibrium", 5, "paper")["None",]$Mean*10^-9
# # #     mirna_new <- paste(unlist(strsplit(mirna,"-")),collapse=".")
# # #     del_g <- dG.df[canonical.sites,mirna_new]
# # #     kds_A <- kds[c(1,2)]
# # #     kds_noA <- kds[c(3, 4)]
# # #     total_x <- c(total_x, kds_A)
# # #     total_y <- c(total_y, kds_noA)
# # #     # b <- optimize(GetConstant,c(-10,20),kds=kds_rel, del_g=del_g)$minimum
# # #     # print(b)
# # #     # y <- k.R*k.T*log(x) + b
# # #     if (mirna == "miR-1") {
# # #       plot(kds_A*1e9, kds_noA*1e9, col=kMirnaColors[mirna],log='xy', xlim = c(1e-3, 5e-1), ylim = c(1e-3, 5e-1), pch = c(1, 2))    
# # #     } else {
# # #       points(kds_A*1e9, kds_noA*1e9, col=kMirnaColors[mirna], pch = c(1, 2))    

# # #     }
# # #     x <- 10^(seq(-4,1, length = 10))
# # #     print(mean(kds_noA/kds_A))
# # #     y <- mean(kds_noA/kds_A)*x
# # #       lines(x, y, lty = 2, col = kMirnaColors[mirna])

# # #     # title(main=mirna)
# # #     # lines(x,y)
# # #     # plot(kds_rel, del_g,col=kSiteColors[canonical.sites,],ylim = c(-15, 5), log='x')
# # #   # lines(x,y)
# # # }
# # #       lines(x, x, lty = 2, col = "gray")

# #   text(0.5, -10, round(cor(total_x, total_y)^2, 2))

# #   legend("bottomright", legend = mirnas.all, col=kMirnaColors[mirnas.all], bty="n",pch=19)


# # }

# # plotkdsvsEnergy()
# # dev.copy2pdf(file="2017_Paper/Figure2E_raw_v6.pdf",useDingbats=FALSE)
# # break
# # MakeFigure4("miR-1","7mer-m8")

# # data.I <- read.table(file = "/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/structures_bp_prob/8mer/I_0-0.txt",
# #                      sep = "\t", header=FALSE, skip = 1)
# # data.A_save <- as.data.frame(
# #   fread(file = "/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/structures_bp_prob/8mer/4_0-0.txt",
# #                      sep = "\t", 
# #                      header=FALSE,
# #                      skip = 1,
# #                      stringsAsFactors=FALSE),
# #   colClasses=c("integer", rep("character", 4), rep("float", 84)))

# # reads.I <- read.table(file = "/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/reads_by_site/8mer/I_0-0.txt",skip=1,stringsAsFactors=FALSE)


# # names(data.I) <- read.table(file = "/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/structures_bp_prob/8mer/I_0-0.txt",
# #                      sep = "\t", header=FALSE, nrow=1)[1:89]
# # names(data.A) <- names(data.I)

# tSNE <- function(data_probs) {
#   len <- nrow(data_probs)
#   # plot(1, type = "n", xlim = c(1, 38), ylim = c(0, 1))
#   output <- matrix(NaN,nrow=len,ncol = 38)
#   for (row in seq(len)) {
#     mir_index <- as.numeric(data_probs[row,1] + 1)
#     range_initial <- c(mir_index - 15, mir_index + 15 + 7)
#     seed_positions <- data_probs[row, (range_initial[1]+5):(range_initial[2]+5)]
#     output[row,1:38] <- as.numeric(unlist(seed_positions))
#   }
#   return(output)
# }

# # data.A <- data.A_save[sample(1:nrow(data.A_save),10),]
# # inds.I <- sample(1:nrow(data.I),5000)
# # inds.A <- sample(1:nrow(data.A),5000)
# # window.I <- tSNE(data.I[inds.I,])
# # cols.I <- apply(data.I[inds.I,],1, function(row){
# #   return(GetColorFunction(paste0(row[2:5],collapse="")))
# #   })
# # new.I <- cbind(window.I, rep(rgb(0,0,1, alpha = 0.5), nrow(window.I)))
# # window.A <- tSNE(data.A[inds.A,])
# # cols.A <- apply(data.I[inds.A,],1, function(row){
# #   return(GetColorFunction(paste0(row[2:5],collapse="")))
# #   })
# # new.A <- cbind(window.A, rep(rgb(1,0,0, alpha = 0.5), nrow(window.A)))
# # merge <- data.frame(data = rbind(window.I, window.A),
# #   colors = c(cols.I, cols.A), type=c(rep("I", nrow(window.I)), rep("A", nrow(window.A))),
# #   stringsAsFactors=FALSE)
# # merge <- merge[!duplicated(merge[,1:38]),]
# # print(3159)
# # set.seed(42)
# # tsne.data <- Rtsne(merge[,1:38],verbose=TRUE)
# # # tsne.A <- Rtsne(unique(window.A))
# # # print(3162)
# # # colors.I <- sapply(unique(window.I), function(row) {
# # #   return(GetColorFunction(paste0(row[2:5], collapse="")))
# # #   })

# # # colors.A <- sapply(unique(window.A), function(row) {
# # #   return(GetColorFunction(paste0(row[2:5], collapse="")))
# # #   })
# # # merged.colors <- merge[rownames(unique(merge[,-ncol(merge),])),]$colors
# # I.inds <- which(merge$type == "I")
# # A.inds <- which(merge$type == "A")

# # # dev.new(xpos = 20, ypos = 20, width = 10, height = 5)
# # # par(mfrow = c(1, 2))


# # plot(tsne.data$Y[I.inds,], type = "p",col=merge$colors[I.inds],pch=1,lwd=1)
# # plot(tsne.data$Y[A.inds,], type = "p",col=merge$colors[A.inds],pch=1,lwd=1)

# # plot(tsne.A$Y)
# # graphics.off()
# # # MakeFigure1()
# # # MakeFigure2()

# # # MakeSupplementalFigure2()

# # MakeFigure3()




# GetAllFlanks <- function(start, stop, sitelist, log.residual=FALSE, combined.input=TRUE) {
#   output <- GetAllMirnaFlankKds("let-7a", "equilibrium", start, stop, sitelist, log.residual=log.residual, combined.input=combined.input)
#   colnames(output) <-sapply(colnames(output), function(name) {paste0("let-7a_", name)})
#   for (mirna in c("miR-1", "miR-155", "miR-124", "lsy-6")) {
#     output.temp <- GetAllMirnaFlankKds(mirna, "equilibrium", start, stop, sitelist, log.residual=log.residual, combined.input=combined.input)
#     colnames(output.temp) <-sapply(colnames(output.temp), function(name) {paste0(mirna, "_", name)})

#     output <- cbind(output, output.temp)
#   }
#   output <- log10(output)
#   output <- t(t(output) - colMeans(output, na.rm =TRUE))
#   return(output)

# }

# MakeKdandSequenceTable <- function(mirna, experiment, start, stop, sitelist,
#                                    log.residual=FALSE, combined.input=TRUE,
#                                    nosite=TRUE){
#   kds <- GetKds(mirna, experiment, start, stop, sitelist,
#                 log.residual=log.residual, combined.input=combined.input,
#                 nosite=nosite)
#   print(kds)
#   data <- GetSitesXCountsCombined(mirna, experiment, start, stop, sitelist)
#   print(data[,1:2])
#   out <- data.frame(seq = data[,1], kd = as.numeric(kds[rownames(data)]))
#   rownames(out) <- rownames(data)
#   return(out)
# }

# MakeSiteKdBeeswarms <- function(mirna, experiment, start, stop, sitelist,
#                                 site_list=NULL,flankdata=NULL,
#                                 colorByPoints=FALSE,
#                                 sitelist.print=FALSE, log.residual=FALSE,
#                                 combined.input=TRUE, nosite=TRUE) {

#   dev.new(xpos = 20, ypos = 20, height = 8, width = 12)
#   if (length(site_list) == 0) {
#     site_list <- rownames(data)
#   } else if (class(site_list) == "character") {
#     site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
#                                             "computation/AgoRBNS/",
#                                             "AssignSiteTypes/sites.", mirna,
#                                             "_", site_list, ".txt"),
#                               stringsAsFactors=FALSE)[,1], "None")
#   } else {
#     site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
#   }
#   num.sites <- length(site_list)
#   # Get kds for all site-types of the mirna.
#   site.kds <- GetKds(mirna, experiment, start, stop, sitelist,
#                      log.residual=log.residual, combined.input=combined.input,
#                      nosite=nosite)
#   print(site.kds)
#   site.kds <- site.kds[1 : (length(site.kds) - 2)]  
#   if (nosite == FALSE) {
#     site.kds <- c(site.kds, 1)
#     names(site.kds)[length(site.kds)] <- "None"
#   }
#   # # Get flanking kds for the mirna.
#   if (length(flankdata) == 0){
#       flank.kds <- GetAllMirnaFlankKds(mirna, experiment, start, stop, sitelist,
#                                    log.residual=log.residual,
#                                    combined.input=combined.input,nosite=nosite)
#   } else {
#     flank.kds <- flankdata
#   }
#   site.kds <- site.kds[which(colSums(is.na(flank.kds))!=256)]
#   site.kds <- site.kds[site_list[-length(site_list)]]
#   site.kds <<- site.kds
#   site_list <<- site_list
#   flank.kds.trim <- flank.kds[,site_list[-length(site_list)]]
#   flank.kds.trim <- flank.kds.trim[,order(site.kds)]
#   flank.kds.trim <<- flank.kds.trim
#   kds.flanks <- c(flank.kds.trim)
#   data.sites <- rep(colnames(flank.kds.trim), each=nrow(flank.kds.trim))
#   data.ranks <- rep(num.sites-seq(num.sites-1), each=nrow(flank.kds.trim))
#   data.colors <- rep(sapply(rownames(flank.kds), GetColorFunction, alpha=0.7),
#                      ncol(flank.kds.trim))
#   flanks.df <- data.frame(kds=log10(kds.flanks), rank=data.ranks, sites=data.sites,
#                           cols=data.colors,
#                            stringsAsFactors=FALSE)
#   print(unique(flanks.df$rank))
#   print(unique(flanks.df$site))

#   print(unique(flanks.df$rank))
#   print(unique(flanks.df$site))
#   flanks.df <<- flanks.df
#   par(kPlotParameters)
#   boxplot(kds ~ rank,
#           data       = flanks.df,
#           axes       = FALSE,
#           horizontal = TRUE,
#           outline    = FALSE,
#           xlim       = c(0, num.sites+5),
#           ylim       = rev(c(-4,1)))
#   print(c(0,num.sites+5))
#   title(main = mirna,
#         line = -2,
#         adj  = 0.1)
#   title(xlab = expression(K[D]))
#   if (colorByPoints == TRUE){
#     beeswarm(kds ~ rank,
#              data       = flanks.df,
#              add        = TRUE,
#              method     = "swarm",
#              corral     = "random",
#              pch        = 1,
#              horizontal = TRUE,  
#              pwcol = cols)
#   } else {
#     beeswarm(kds ~ rank,
#          data       = flanks.df,
#          add        = TRUE,
#          method     = "swarm",
#          corral     = "random",
#          pch        = 20,
#          horizontal = TRUE,  
#          col        = kSiteColors[rev(unique(flanks.df$sites)), ])
#   }
#     # points(log10(site.kds[order(site.kds)]),num.sites - seq(num.sites-1),col=kSiteColors[names(site.kds),][order(site.kds)],cex=3)

#   ymin=0.0001
#   ymax=10
#   ys <- c(sapply(seq(floor(min(flanks.df$kds, na.rm=TRUE)),
#                      ceiling(max(flanks.df$kds, na.rm=TRUE))), function(x) seq(10)*10^x))
#   ys <- ys[ys >= ymin & ys <= ymax]
#   ymin <- min(ys)
#   ymax <- max(ys)
#   yl <- seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#   axis(1, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=0, lwd=0, cex.axis=1.5, padj=0.2)
#   axis(1, at=log10(ys), labels=FALSE,
#        pos=0, lwd=2)
#   points(rep(-3.4, num.sites-1),num.sites - seq(num.sites-1),col= kSiteColors[names(site.kds),][order(site.kds)],cex=1.5,pch=19)

#   for (ind in seq(num.sites)) {
#     site = unique(flanks.df$sites)[ind]
#     text(-3.5,num.sites - ind, labels=site, adj=0, cex=1.2, col= "black")
#   }
# }

# MakeSiteBarPlotsOrig <- function(mirna, experiment, start, stop, sitelist, site_list, 
#                                 colorByPoints=FALSE,
#                                 sitelist.print=FALSE, log.residual=FALSE,
#                                 combined.input=TRUE) {
#   # Get kds for all site-types of the mirna.
#   dev.new(xpos = 20, ypos = 20, height = 8, width = 12)
#   site.kds <- GetKdsErrorOrig(mirna, experiment, start, stop, sitelist)
#   print(site.kds)
#       site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", site_list,".txt"),
#   stringsAsFactors=FALSE)[,1], "None")
#   print(site.kds[1:10,])
#   print(site_list)
#   print(rownames(site.kds))
#   print(site.kds)
#   site.kds <- site.kds[which(rownames(site.kds) %in% site_list),]
#   site.kds <- site.kds / c(site.kds["None",])
#   flanks.df <- data.frame(kds=site.kds[,2],upper=site.kds[,2],lower=site.kds[,3], site = rownames(site.kds),stringsAsFactors=FALSE)
#   flanks.df <- flanks.df[order(flanks.df$kds),]
#   print(flanks.df)
#   num.sites <- nrow(site.kds)
#   print(dim(flanks.df))
#   par(kPlotParameters)
#   xs <- site.kds
#   ys <- nrow(site.kds) - seq(nrow(site.kds)) + 1
#   xs <<- xs
#   ys <<- ys
#   print(length(flanks.df$kds))
#   print(kSiteColors[flanks.df$site,])
#   print(length(nrow(site.kds) - seq(nrow(site.kds)) + 1))
#   plot(flanks.df$kds,nrow(site.kds) - seq(nrow(site.kds)) + 1,
#           col = "white",
#           axes    = FALSE,
#           log = 'x',
#           pch = 20,
#           cex = 2,
#           ylim       = c(0, num.sites+5),
#           xlim       = rev(c(0.00001, 1.1)))
#  arrows(flanks.df$upper, nrow(site.kds) - seq(nrow(site.kds)) + 1,
#         flanks.df$lower, nrow(site.kds) - seq(nrow(site.kds)) + 1, length=0.07, angle=90, code=3)

#   title(main = mirna,
#         line = -2,
#         adj  = 0.1)
#   title(xlab = expression(K[D]))
#  points(flanks.df$kds,nrow(site.kds) - seq(nrow(site.kds)) + 1,
#           bg = kSiteColors[flanks.df$site,],
#           col = "black",
#           pch = 21,
#           cex = 1.5)
#   ymin=0.0001
#   ymax=1
#   ys <- c(sapply(seq(log10(ymin),
#                      log10(ymax)), function(x) seq(10)*10^x))
#   # ys <- ys[ys >= ymin & ys <= ymax]
#   ymin <- min(ys)
#   ymax <- max(ys)
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
#   axis(1, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=0, lwd=0, cex.axis=1.5, padj=0.2)
#   axis(1, at=ys, labels=FALSE,
#        pos=0, lwd=2)
# print(unique(flanks.df$rank))
# print(unique(flanks.df$sites))
# print(flanks.df[1:10,])
# points(rep(0.00006, nrow(site.kds)),nrow(site.kds) - seq(nrow(site.kds)) + 1,col="black",bg=kSiteColors[flanks.df$site,],cex=1.5,pch=21)
# sapply(seq(num.sites), function(ind) {
#   site = flanks.df$site[ind]

#   text(0.00005,num.sites - ind + 1, labels=site, adj=0, cex=1.2, col= "black")
#   })
# }






# PlotFlankLinearModel <- function(flanks, site_fit, site_y,  int2=FALSE, int3=FALSE) {

#   flanks <- log10(flanks)
#   flanks <- t(t(flanks) - colMeans(flanks, na.rm=TRUE))
#   fit <- LinearModelInputFlanks(flanks, which(colnames(flanks) == site_fit))
#   if (int3 == TRUE) {
#     fit2 <- lm(I ~ f5p.i*f5p.o*f3p.i*f3p.o, data=fit)      
#   } else if (int2 == TRUE) {
#     fit2 <- lm(I ~ f5p.i*f5p.o + f3p.i*f3p.o, data=fit)      
#   } else {
#     fit2 <- lm(I ~ f5p.i + f5p.o + f3p.i + f3p.o, data=fit)  
#   }
#   plot(  fit2$fitted.values,  flanks[,site_y], xlab = site_fit, pch = 19, cex = 2.0, ylab = site_y, col = kNucleotideColors[fit$f5p.o])
#   points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 1.5, col = kNucleotideColors[fit$f5p.i])
#   points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 1.0, col = kNucleotideColors[fit$f3p.i])
#   points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 0.5, col = kNucleotideColors[fit$f3p.o])
#   cor_text <- round(
#                 cor(fit2$fitted.values, flanks[, site_y], use="complete"),
#                 digits = 3
#               )

#   title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)

# }

# PlotFlankLinearModelAverage <- function(flanks, site_y, int2=FALSE, int3=FALSE) {

#   flanks <- log10(flanks)
#   flanks <- t(t(flanks) - colMeans(flanks, na.rm=TRUE))
#   flanks <- cbind(rowMeans(flanks, na.rm=TRUE), flanks)
#   fit <- LinearModelInputFlanks(flanks, 1)
#   if (int3 == TRUE) {
#     fit2 <- lm(I ~ f5p.i*f5p.o*f3p.i*f3p.o, data=fit)      
#   } else if (int2 == TRUE) {
#     fit2 <- lm(I ~ f5p.i*f5p.o + f3p.i*f3p.o, data=fit)      
#   } else {
#     fit2 <- lm(I ~ f5p.i + f5p.o + f3p.i + f3p.o, data=fit)  
#   }
#   # print(fit)
#   # print(fit2)
#   # # print(flanks[, site_y])
#   plot(  fit2$fitted.values,  flanks[,site_y], xlab = "average", pch = 19, cex = 2.0, ylab = site_y, col = kNucleotideColors[fit$f5p.o])
#   points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 1.5, col = kNucleotideColors[fit$f5p.i])
#   points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 1.0, col = kNucleotideColors[fit$f3p.i])
#   points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 0.5, col = kNucleotideColors[fit$f3p.o])
#   cor_text <- round(
#                 cor(fit2$fitted.values, flanks[, site_y], use="complete"),
#                 digits = 3
#               )

#   title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)

# }






# MakeWithWithoutProteinScatter <- function(mirna, experiment, start, stop, sitelist) {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
#     colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#       colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#         return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#       }
#     )
#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       k.c.stockago.1 <- 10^params["AGO"]
#       bgs.1 <- 10^params["bg"]

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.2 <- stockago[mirna,experiment]

#       bgs.2 <- 10^params["bg"]

#     print(kds.1)
#     print(kds.2)
#     print(k.c.stockago.1)
#     print(k.c.stockago.2)
#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
#     print(xmin)
#     print(xmax)
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(c(kds.1, bgs.1, k.c.stockago.1),
#       c(kds.2, bgs.2, k.c.stockago.2),
#       col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
#       pch = 19,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#         segments(xmin, xmin, xmax, xmax, lty = 2)

#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = "Fixed protein parameters (nM)", cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = "Fit protein parameters (nM)", cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
# }


# MakeScatterAcrossWindows <- function(mirna, experiment, lim1, lim2, sitelist, method) {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, lim1, lim1, sitelist)
#     colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#       colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#         return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#       }
#     )
#     if (method == "free_protein") {

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       k.c.stockago.1 <- 10^params["AGO"]
#       bgs.1 <- 10^params["bg"]

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       k.c.stockago.2 <- 10^params["AGO"]
#       bgs.2 <- 10^params["bg"]

#       } else {
#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.1 <- stockago[mirna,experiment]

#       bgs.1 <- 10^params["bg"]


#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.2 <- stockago[mirna,experiment]

#       bgs.2 <- 10^params["bg"]
# }
#     print(kds.1)
#     print(kds.2)
#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
#     print(xmin)
#     print(xmax)
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(c(kds.1, bgs.1, k.c.stockago.1),
#       c(kds.2, bgs.2, k.c.stockago.2),
#       col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
#       pch = 19,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#     segments(xmin, xmin, xmax, xmax, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = paste0(lim1, " parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0(lim2, " parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
# }

# MakeScatterAcrossOrder <- function(mirna, experiment, lim, sitelist1, sitelist2, method) {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, lim, lim, sitelist1)
#     colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#       colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#         return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#       }
#     )
#     sitesXcounts1 <- sitesXcounts
#         sitesXcounts <- GetSitesXCounts(mirna, experiment, lim, lim, sitelist2)
#     colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#       colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#         return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#       }
#     )
#     sitesXcounts2 <- sitesXcounts
#     print(sitesXcounts2)
#     if (method == "free_protein") {

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
#       k.c.stockago.1 <- 10^params["AGO"]
#       bgs.1 <- 10^params["bg"]

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
#       k.c.stockago.2 <- 10^params["AGO"]
#       bgs.2 <- 10^params["bg"]

#       } else {
#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.1 <- stockago[mirna,experiment]

#       bgs.1 <- 10^params["bg"]


#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.2 <- stockago[mirna,experiment]

#       bgs.2 <- 10^params["bg"]
# }
#     names(kds.1) <- rownames(sitesXcounts1)
#     names(kds.2) <- rownames(sitesXcounts2)
#     kds.2 <- kds.2[names(kds.1)]
#     print(kds.1)
#     print(kds.2)
#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
#     print(xmin)
#     print(xmax)
#     print(bgs.1)
#     print(bgs.2)
#     print(k.c.stockago.1)
#     print(k.c.stockago.2)
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(c(kds.1, bgs.1, k.c.stockago.1),
#       c(kds.2, bgs.2, k.c.stockago.2),
#       col = kSiteColors[c(rownames(sitesXcounts1), "bg", "Ago"),],
#       pch = 19,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#     segments(xmin, xmin, xmax, xmax, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = paste0(sitelist1, " parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0(sitelist2, " parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts1), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
# }






# MakeSiteIterationPlotKinetics <- function(out,method,mirna, num.kds, num.bgs, colors = FALSE) {
#   setEPS()
#   postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/koffs/",
#     mirna,"/", method, ".eps"))
#   par(kPlotParameters)
#   x = seq(1,dim(out)[1],length = max(ncol(out),1000))
#   out <- out[x,]
#   # print(out)
#   probs = out[ ,"-logp"]
#   out.print <- out[ , seq(dim(out)[2] - 1)]
#   out.print.koffs <- Logistic(out[,1:num.kds],max = 10)
#   out.print.percentages <- Logistic(out[,(num.kds+1):(2*num.kds)],max = 1)
#   out.print.setpoint <- 10^out[,2*num.kds + 1]
#   out.print.bgs <- 10^out.print[,(2 * num.kds + 2) : (2 * num.kds + num.bgs + 1)]

#   out.print <- cbind(out.print.koffs,out.print.percentages,out.print.setpoint, out.print.bgs)
#   ys <- 10^c(max(floor(log10(min(out.print))), -6), ceiling(log10(max(out.print))))
#   probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])


#   plot(x    = x,
#        y    = 10^probs.scaled,
#        log  = "y",
#        axes = FALSE,
#        ann = FALSE, 
#        type = "l",
#        col = "white",
#        ylim = ys)
#   axis(side = 1,
#        at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
#        pos  = ys[1],
#        tck  = -0.01,
#        labels = FALSE)
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
#   axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
#        labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=1, lwd=2)
#   title(main = mirna, font=1)
#   axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
#        pos=1)
#   if (colors != FALSE) {

#   cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
#   names(cols) <- colnames(out.print)
#   cols["AGO"] <- "grey"

#   sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=cols[name])
#     })
#   } else  {
#     sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=kSiteColors[name,])
#     })
#   }
#   lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
#   dev.off()
# }

# centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#                     "11mer-m4.14", "12mer-m4.15")


# kSiteColors <- read.table(
#   "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_info.txt",
#   row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 1, drop = FALSE]





# Norm <- function(vector) {
#   return(vector/ sum(vector))
# }

# Cumul <- function(vector) {
#   norm <- Norm(vector)
#   tot <- 0
#   out <- sapply(norm, function(x){
#     tot <<- tot + x
#     return(tot)

#     })
#   return(out)
# }

# Logistic <- function(vector, max) {
#   return(max / (1 + exp(-vector)))
# }

# Logit <- function(vector, max) {
#   return(-log(max / vector - 1))
# }


# MakeIterationPlot <- function(out,type,extension = "") {
#   setEPS()
#   postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                          "/figures/kds/", mirna, "/", type,
#                          "/iterations/", site, "_", k.c.stockago, extension,".eps"))
#   par(kPlotParameters)
#   x = seq(dim(out)[1])
#   probs = out[ ,"-logp"]
#   out.print <- out[ , seq(dim(out)[2] - 1)]
#   ys <- 10^c(floor(min(out.print)), ceiling(max(out.print)))
#   probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])
#   out.print <- 10^out.print
#   plot(x , 10^probs.scaled, log='y', axes=FALSE, type="l", ylim=ys,
#        lwd=2, ann=FALSE,
#        col="black")
#   title(main = mirna, line=-1, adj=0.1)
#   title(main = site, col.main=kSiteColors[site,], line=-2.5, adj=0.1)

#   title(xlab = "Iteration")
#   title(ylab = "Parameter values (nM)")
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))), pos=ys[1], lwd=2,
#        labels=FALSE, tck=-0.01)
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
#   axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
#        labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=1, lwd=2, hadj=0.8)
#   axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
#        pos=1, lwd=2)
#   sapply(colnames(out.print), function(name) {
#           lines(x, out.print[, name], lwd=2, col=GetColorFunction(name))
#         }
#         )
#   lines(x, 10^probs.scaled, type="l", col="black")
#   dev.off()
# }



# WriteIterationFile <- function(out,extension="") {
#   out.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                      "/equilibrium/kds/", site,"_flanking_",
#                      k.c.stockago, extension,".txt")
#   write.table(file=out.file, out, sep="\t", quote=FALSE, row.names=FALSE,
#                 col.names=TRUE)
# }

# WriteFinalParameterFile <- function(out,extension="") {
#   out.final <- out[dim(out)[1], ]
#   names(out.final) <- colnames(out)
#   out.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                            "/equilibrium/kds/final_", site,
#                            "_flanking_", k.c.stockago,extension, ".txt")
#   write.table(file=out.file, out.final, sep="\t", quote=FALSE, row.names=TRUE,
#               col.names=FALSE)
# }












# # FIXED FOR PAPER
# PlotSiteKdOptimization <- function(out, specific_string, mirna, num.kds,
#                                    num.bgs, colors = FALSE) {
#   # Assign file name to the figure.
#   setEPS()
#   filename <- paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/",
#                          "2017_Paper/Kd_fits/", mirna, "_", specific_string, ".eps")
#   postscript(file=filename)
#   # Load standard plot settings.
#   par(kPlotParameters)
#   # Define x as either the list of all rows, or alternatively 1000 equally
#   # spaced rows.
#   x <- seq(1, nrow(out), length=min(nrow(out), 1000))
#   out <- out[x, ]
#   # Split output matrix into loglikelihood and parameters:
#   probs <- out[, "-logp"]
#   pars <- out[, seq(ncol(out) - 1)]
#   kds <- Logistic(pars[, 1:num.kds], max = 10)
#   bgs <- 10^pars[,(num.kds + 1) : (num.kds + num.bgs)]
#   AGO <- 10^pars[,(num.kds + num.bgs + 1) : ncol(pars)]
#   out.print <- cbind(kds, bgs, AGO)
#   ys <- c(10^-10, 10^5)
#   probs.scaled <- probs / probs[1]

#   plot(x    = x,
#        y    = 10^probs.scaled,
#        log  = "y",
#        axes = FALSE,
#        type = "l",
#        col = "white",
#        ylim = ys)
#   axis(side = 1,
#        at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
#        pos  = ys[1], labels = FALSE,
#        tck  = -0.01)
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
#   axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
#        labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=1, lwd=2)
#   title(main = mirna, font=1)
#   axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
#        pos=1)
#   if (colors != FALSE) {

#   cols <- colors
#   names(cols) <- colnames(out.print)
#   cols["AGO"] <- "grey"

#   sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=cols[name])
#     })
#   } else  {
#     sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=kSiteColors[name,])
#     })
#   }
#   lines(x, ys[1]*(ys[2]/ys[1])^probs.scaled, type="l", lwd=2, col="black")
#   dev.off()
# }


# PlotSiteFlankKdOptimization <- function(out, specific_string, mirna, num.kds,
#                                    num.bgs, colors) {
#   # Assign file name to the figure.
#   setEPS()
#   postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/",
#                          "2017_Paper/Kd_fits/", mirna, "_", specific_string, ".eps"))
#   # Load standard plot settings.
#   par(kPlotParameters)
#   # Define x as either the list of all rows, or alternatively 1000 equally
#   # spaced rows.
#   x <- seq(1, nrow(out), length=min(nrow(out), 1000))
#   out <- out[x, ]

#   # Split output matrix into loglikelihood and parameters:
#   probs <- out[, "-logp"]
#   pars <- out[, seq(ncol(out) - 1)]
#   out.print <- Logistic(pars, max = 10)

#   ys <- 10^c(max(floor(log10(min(out.print))), -10), ceiling(log10(max(out.print))))
#   probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])

#   plot(x    = x,
#        y    = 10^probs.scaled,
#        log  = "y",
#        axes = FALSE,
#        type = "l",
#        col = "white",
#        ylim = ys)
#   axis(side = 1,
#        at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
#        pos  = ys[1], labels = FALSE,
#        tck  = -0.01)
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
#   axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
#        labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=1, lwd=2)
#   title(main = mirna, font=1)
#   axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
#        pos=1)

#   cols <- colors
#   print(dim(out.print))
#   print(length(cols))
#   names(cols) <- colnames(out.print)

#   sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=cols[name])
#     })
#   lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
#   dev.off()
# }







# MakeEquilibriumEnrichmentPlots <- function(mirna, experiment, start, stop, sitelist, site_list=NULL, combined.input=TRUE, log.residual=FALSE, bgoff=FALSE, model=TRUE, connected.points=FALSE,nosite=FALSE) {
#   dev.new(xpos = 20, ypos = 20, height = 6.9, width = 8.157)
#   params <- GetKds(mirna, experiment, start, stop, sitelist, combined.input=combined.input, log.residual=log.residual, scaled=FALSE,nosite=nosite)
#   if (combined.input == TRUE) {
#     sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop, sitelist)
#   } else {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
#   }
#   sitesXcounts <- sitesXcounts[,-1]

#   # Remove any rows with no reads whoatsoever.
#   sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
#   # Remove any rows for which there are no input reads:
#   sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
#   if (nosite == FALSE){
#     kds <- c(unlist(Logistic(params[1:(nrow(sitesXcounts)-1)], 10)),1)
#     } else{
#       kds <- Logistic(params[1:nrow(sitesXcounts)], 10)
#     }
#   names(kds) <- rownames(sitesXcounts)

#   bgs <- rep(10^params["bg"], 5)
#   k.c.stockago <- 10^params["AGO"]

#   c.I.tots <- Norm(sitesXcounts[,1])*100
#   names(c.I.tots) <- rownames(sitesXcounts)
#   data <- sitesXcounts[,2:6]
#   names(c.I.tots) <- rownames(data)
#   x_model <- seq(min(as.numeric(colnames(data))) * 0.7, max(as.numeric(colnames(data))) / 0.7,
#                  length=100)

#   c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
#   colnames(c.totals) <- colnames(x_model)
#   rownames(c.totals) <- rownames(data)

#   c.agos <- sapply(x_model, function(x) {
#       as.numeric(x) * k.c.stockago / 100
#     }
#   )

#   c.bounds <- as.matrix(
#     sapply(c.agos, function(x) {
#         return(GetBoundRNA(kds, c.I.tots, x))
#       }
#     )
#   )
#   print("percent free")
#   print((c.agos - colSums(c.bounds) )/ c.agos)
#   c.frees <- c.totals - c.bounds
#   c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))

#   if (bgoff == TRUE) {
#     c.all <- c.bounds
#   } else {
#     c.all <- c.bounds + c.bgs
#   }

#   c.final <- data.frame(t(t(c.all) / colSums(c.all)))
#   rownames(c.final) <- rownames(data)
#   x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000/2
#   print(x)
#   y <- c(1,1,1,1,1)
#   print(x)
#   sites.norm <- Norm(c.I.tots)
#   data.norm <- t(t(data)/colSums(data))

#   data.R <- data.norm/(sites.norm)
#   model.R <- c.final/(sites.norm)

#   xmin <- min(x)*0.3
#   xmax <- max(x)*3
#   ymin <- 0.2
#   ymax <- 300
#   yextension <- (ymax/ymin)
#   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   xs <- xs[xs >= xmin & xs <= xmax]
#   xmin <- min(xs)
#   xmax <- max(xs)
#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))
#   ys <- ys[ys >= ymin & ys <= ymax]
#   ymin <- min(ys)
#   ymax <- max(ys)
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))

#   plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.5), ylim=c(ymin, ymax), type="l",
#      col="white", axes=FALSE, ann=FALSE)        
#   # Generate tickmarks for axis.

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin, lwd=2)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#   axis(2, at=ys, labels=FALSE,
#        pos=xmin, lwd=2)

#   title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#   title(xlab = "[AGO2-miRNA] (pM)", cex.lab=1.5, line=2, adj=0.3)
#   title(ylab = "Enrichment", cex.lab=1.5, line=2)

#   if (length(site_list) == 0) {
#     site_list <- rownames(data)
#   } else if (class(site_list) == "character") {
#     site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
#                                             "computation/AgoRBNS/",
#                                             "AssignSiteTypes/sites.", mirna,
#                                             "_", site_list, ".txt"),
#                               stringsAsFactors=FALSE)[,1], "None")
#   } else {
#     site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
#   }

#   legend.names <- rownames(data)[order(kds)]
#   print(order(kds))
#   legend.names <- legend.names[which(legend.names %in% site_list)]
#   ordered_list <- legend.names


#   centered.sites <- c("11mer-m3.13", "12mer-m3.14", "11mer-m4.14",
#                       "12mer-m4.15", "13mer-m3.15", "11mer-m5.15")
#   centered_rename <- c("Centered-m3.13", "Centered-m3.14", "Centered-m4.14",
#                        "Centered.m4-15", "Centered-m3.15", "Centered-m5.15")
#   names(centered_rename) <- centered.sites

#   legend.names <- sapply(legend.names, function(site) {
#                            if (site %in% centered.sites) {
#                              return(centered_rename[site])
#                            } else {
#                              return(site)
#                           }
#                          }
#                          )

#   legend(x=xmax, y=ymax, legend=legend.names, pch=19,
#          col=kSiteColors[ordered_list], cex=1, bty="n", ncol=1)
#   for (name in site_list) {
#     if (connected.points == TRUE) {
#       type = "o"
#     } else {
#       type = "p"
#     }
#     points(x, data.R[name, ], col=kSiteColors[name], type=type, pch=19,
#            cex=1.5, lwd=3)
#     if (model == TRUE) {
#     lines(x_model*k.c.stockago/100*1000/2, model.R[name, ], col=kSiteColors[name], lwd=2)      
#     }
#   }
# }

# MakeOccupancyPlots <- function(mirna, experiment, start, stop, sitelist, site_list=NULL, combined.input=TRUE, log.residual=FALSE, bgoff=FALSE, model=TRUE, connected.points=FALSE,nosite=FALSE) {
#   params <- GetKds(mirna, experiment, start, stop, sitelist, combined.input=combined.input, log.residual=log.residual, scaled=FALSE, nosite=nosite)
#   if (combined.input == TRUE) {
#     sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop, sitelist)
#   } else {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
#   }
#   sitesXcounts <- sitesXcounts[,-1]

#   # Remove any rows with no reads whoatsoever.
#   sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
#   # Remove any rows for which there are no input reads:
#   sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

#   if (nosite == FALSE){
#     kds <- c(unlist(Logistic(params[1:(nrow(sitesXcounts)-1)], 10)),1)
#     } else{
#       kds <- Logistic(params[1:nrow(sitesXcounts)], 10)
#     }
#   names(kds) <- rownames(sitesXcounts)

#   bgs <- rep(10^params["bg"], 5)
#   k.c.stockago <- 10^params["AGO"]
#   print(k.c.stockago)
#   c.I.tots <- Norm(sitesXcounts[,1])*100
#   names(c.I.tots) <- rownames(sitesXcounts)
#   data <- sitesXcounts[,2:6]
#   names(c.I.tots) <- rownames(data)
#   x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100
#   print(x)
#   x_model <- seq(min(as.numeric(colnames(data))) * 0.7, max(as.numeric(colnames(data))) / 0.7,
#                  length=100)

#   c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
#   colnames(c.totals) <- colnames(x_model)
#   rownames(c.totals) <- rownames(data)

#   c.agos <- sapply(x_model, function(x) {
#       as.numeric(x) * k.c.stockago / 100
#     }
#   )
#   c.bounds <- as.matrix(
#     sapply(c.agos, function(x) {
#         return(GetBoundRNA(kds, c.I.tots, x))
#       }
#     )
#   )
#   c.bounds.points <- as.matrix(sapply(x,function(i) {
#     return(GetBoundRNA(kds,c.I.tots,i))
#     }
#     )
#   )
#   c.frees <- c.totals - c.bounds


#   c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))

#   if (bgoff == TRUE) {
#     c.all <- c.bounds
#   } else {
#     c.all <- c.bounds + c.bgs
#   }

#   c.final <- data.frame(t(t(c.all) / colSums(c.all)))
#   rownames(c.final) <- rownames(data)
#   print(x)
#   y <- c(1,1,1,1,1)
#   print(x)
#   sites.norm <- Norm(c.I.tots)
#   data.norm <- t(t(data)/colSums(data))

#   data.R <- data.norm/(sites.norm)
#   model.ago.occupancy <- t(t(c.bounds) / colSums(c.bounds))
#   model.ago.occupancy.points <- t(t(c.bounds.points) / colSums(c.bounds.points))
#   xmin <- min(x)*0.3
#   xmax <- max(x)*3
#   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   xs <- xs[xs >= xmin & xs <= xmax]
#   xmin <- min(xs)
#   xmax <- max(xs)
#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   ymin <- 0
#   ymax <- 0.7
#   ys <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
#   yl <- ys

#   plot(x, y,log='x', xlim=c(xmin, xmax*(xmax/xmin)^0.55), ylim=c(ymin, ymax), type="l",
#      col="white", axes=FALSE, ann=FALSE)        
#   # Generate tickmarks for axis.

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin, lwd=2)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=yl,
#        pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#   axis(2, at=ys, labels=FALSE,
#        pos=xmin, lwd=2)

#   title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#   title(xlab = "[AGO2-miRNA] (pM)", cex.lab=1.5, line=2, adj=0.3)
#   title(ylab = "Enrichment", cex.lab=1.5, line=2)
#   print(site_list)
#   if (length(site_list) == 0) {
#     site_list <- rownames(data)
#   } else if (class(site_list) == "character") {
#     site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", site_list,".txt"),
#   stringsAsFactors=FALSE)[,1], "None")
#   } else {
#     site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
#   }
#   print(site_list)
#   legend.names <- rownames(data)[order(kds)]

#   legend.names <- legend.names[which(legend.names %in% site_list)]
#   ordered_list <- legend.names
#   print(ordered_list)
#   centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#                   "11mer-m4.14", "12mer-m4.15", "13mer-m3.15", "11mer-m5.15")
#   centered_rename <- c("Centered-m3.13", "Centered-m3.14", "Centered-m4.14", "Centered.m4-15", "Centered-m3.15", "Centered-m5.15")
#   names(centered_rename) <- centered.sites
#   print(kSiteColors[site_list, ])
#   print(kSiteColors[ordered_list])

#   legend.names <- sapply(legend.names, function(site){
#                         if (site %in% centered.sites) {
#                           return(centered_rename[site])
#                         } else {
#                           return(site)
#                         }
#                        }
#                        )
#   # legend.names <- unique(legend.names)
#   print(legend.names)
#   legend(x=xmax,y=ymax,legend=legend.names, pch=19, col=kSiteColors[ordered_list], cex=1.1, bty="n", ncol = 1)
#   for (name in site_list) {
#     if (connected.points == TRUE) {
#       type = "o"
#     } else {
#       type = "p"
#     }
#     if (model == TRUE) {
#     lines(x_model*k.c.stockago/100*1000/2, model.ago.occupancy[name, ], col=kSiteColors[name], lwd=2)      
#     points(x,model.ago.occupancy.points[name,],col=kSiteColors[name], pch=19,cex=1.5)
#     }
#   }
# }








# MakeFlankingEquilibriumEnrichmentPlots <- function(mirna, experiment, start,
#                                                    stop, site, sitelist, 
#                                                    site_list=NULL,
#                                                    combined.input=TRUE,
#                                                    log.residual=FALSE,
#                                                    bgoff=FALSE, model=TRUE,
#                                                    connected.points=FALSE,
#                                                    nosite=TRUE) {
#   params <- GetKds(mirna, experiment, start, stop, sitelist,
#                    combined.input=combined.input, log.residual=log.residual,
#                    scaled=FALSE, nosite=nosite)
#   params.flanks <- GetFlankKds(mirna, experiment, site, start, stop, sitelist,
#                                combined.input=combined.input,
#                                log.residual=log.residual,scaled=FALSE,nosite=nosite)

#   if (combined.input == TRUE) {
#     sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop,
#                                             sitelist)
#   } else {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
#   }
#   sitesXcounts <- sitesXcounts[,-1]
#   print(params)
#   print(params.flanks)
#   # Remove any rows with no reads whoatsoever.
#   sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
#   # Remove any rows for which there are no input reads:
#   sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
#   if (nosite==TRUE){
#   kds.s <- unlist(Logistic(params[1:nrow(sitesXcounts)], 10))
#   } else {
#   kds.s <- c(unlist(Logistic(params[1:(nrow(sitesXcounts)-1)], 10)), 1)

#   }
#   names(kds.s) <- rownames(sitesXcounts)
#   s.c <- as.numeric(sitesXcounts[site, ])
#   sfXc <- GetSiteFlanksXCounts(mirna, experiment, site, start, stop, sitelist)
#   sfXc <- sfXc[,-1]

#   colnames(sfXc)[1] <- "I"
#   sfXc <- sfXc[rowSums(sfXc[, 2:6]) > 0,]
#   sfXc <- sfXc[which(sfXc[,1] > 0),,drop=FALSE]

#   sfXc <- round(t((t(sfXc) / colSums(sfXc)) * s.c))
#   sfXc[is.na(sfXc)] <- 0
#   sitesXcounts.sites <- sitesXcounts[rownames(sitesXcounts) != site,]
#   sitesXcounts <- rbind(sitesXcounts, sfXc)
#   sitesXcounts <- sitesXcounts[rownames(sitesXcounts) != site, ]
#   print(s.c)

# # Omit the site kd for which the flanking sites are being fit.
#   kds.s <- kds.s[names(kds.s) != site]
#   kds <- c(kds.s, Logistic(params.flanks, max=10))
#   print(kds)
#   bgs <- rep(10^params["bg"], 5)
#   k.c.stockago <- 10^params["AGO"]
#   print(bgs)
#   print(k.c.stockago)
#   c.I.tots <- Norm(sitesXcounts[,1])*100
#   names(c.I.tots) <- rownames(sitesXcounts)
#   data <- sitesXcounts[,2:6]
#   c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 5, byrow=FALSE)
#   colnames(c.totals) <- colnames(data)
#   rownames(c.totals) <- rownames(data)
#   colors_sites <- kSiteColors[rownames(sitesXcounts)[rownames(sitesXcounts) != site],]
#   colors_flanks <- sapply(rownames(sfXc), GetColorFunction)
#   colors_all <- c(rep("grey", length(kds.s)), colors_flanks)
#   names(colors_all) <- rownames(sitesXcounts)
#   x_model <- seq(min(as.numeric(colnames(data))) * 0.7, max(as.numeric(colnames(data))) / 0.7,
#                  length=100)

#   c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
#   colnames(c.totals) <- colnames(x_model)
#   rownames(c.totals) <- rownames(data)

#   c.agos <- sapply(x_model, function(x) {
#       as.numeric(x) * k.c.stockago / 100
#     }
#   )

#   c.bounds <- as.matrix(
#     sapply(c.agos, function(x) {
#         return(GetBoundRNA(kds, c.I.tots, x))
#       }
#     )
#   )
#   c.frees <- c.totals - c.bounds
#   c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))

#   if (bgoff == TRUE) {
#     c.all <- c.bounds
#   } else {
#     c.all <- c.bounds + c.bgs
#   }

#   c.final <- data.frame(t(t(c.all) / colSums(c.all)))
#   rownames(c.final) <- rownames(data)
#   x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000
#   y <- c(1,1,1,1,1)

#   sites.norm <- Norm(c.I.tots)
#   data.norm <- t(t(data)/colSums(data))

#   data.R <- data.norm/(sites.norm)
#   model.R <- c.final/(sites.norm)
#   data.R <<- data.R
#   model.R <<- model.R
#   xmin <- min(x)*0.3
#   xmax <- max(x)*3
#   ymin <- 0.2
#   ymax <- 300
#   yextension <- (ymax/ymin)
#   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   xs <- xs[xs >= xmin & xs <= xmax]
#   xmin <- min(xs)
#   xmax <- max(xs)
#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))
#   ys <- ys[ys >= ymin & ys <= ymax]
#   ymin <- min(ys)
#   ymax <- max(ys)
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))

#   plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.55), ylim=c(ymin, ymax), type="l",
#      col="white", axes=FALSE, ann=FALSE)        
#   # Generate tickmarks for axis.

#   axis(1, at=xl,
#        labels=sapply(xl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin, lwd=2)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=sapply(yl, function(name) {
#          eval(substitute(expression(10^x), list(x=log10(name))))
#        }),
#        pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#   axis(2, at=ys, labels=FALSE,
#        pos=xmin, lwd=2)

#   title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#   title(xlab = "[AGO2-miRNA] (pM)", cex.lab=1.5, line=2, adj=0.3)
#   title(ylab = "Enrichment", cex.lab=1.5, line=2)
#   print(site_list)
#   if (length(site_list) == 0) {
#     site_list_real <- rownames(data)
#   } else if (class(site_list) == "character") {
#     site_list_real <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", site_list,".txt"),
#   stringsAsFactors=FALSE)[,1], "None")
#   } else {
#     site_list_real <- c(rownames(sitesXcounts.sites)[order(kds.s)][1:site_list], "None")
#   }
#   site_list_real <- site_list_real[site_list_real!=site]
#   print(site_list_real)
#   legend.names <- rownames(sfXc)

#   ordered_list <- legend.names
#   print(site_list_real)
#   legend.names <- unique(legend.names)
#   print(legend.names)
#   legend(x=xmax,y=ymax,legend=legend.names, text.font=2, text.col= sapply(legend.names,GetColorFunction), cex=0.9, bty="n", ncol = 8)
#   for (name in c(site_list_real,rownames(sfXc))) {
#     print(name)
#     if (connected.points == TRUE) {
#       type = "o"
#     } else {
#       type = "p"
#     }
#     points(x, data.R[name, ], col=colors_all[name], type = type, pch=19, cex=1.5, lwd=3)
#     if (model == TRUE) {
#     lines(x_model*k.c.stockago/100*1000, model.R[name, ], col=colors_all[name], lwd=2)      
#     }
#   }
# }

# PlotAllCenteredKds <- function(start, stop, sitelist, combined.input=TRUE,
#                                log.residual=FALSE, nosite=TRUE) {
#   # This pulls the last two digits out of the sitelist
#   # ("13" from "centered13", etc.)
#   split.sitelist <- unlist(strsplit(sitelist, split=""))
#   split.sitelist.last <- split.sitelist[(nchar(sitelist)-1):nchar(sitelist)]
#   split.sitelist.last <<- split.sitelist.last
#   length_sitelist <- paste(split.sitelist.last,collapse="")

#   kds_centered <- matrix(unlist(sapply(mirnas.all, function(mirna) {
#     kds <- unlist(GetKds(mirna,"equilibrium", start, stop, sitelist,
#                   combined.input=combined.input, log.residual=log.residual,
#                   nosite=nosite))
#     print(kds)
#     kds <- c(kds["8mer"],kds["6mer"], kds["6mer-m8"],kds[grep("\\.",names(kds)),drop=FALSE],kds["None"])
#     kds <- kds[grep("23", names(kds), invert=TRUE),drop=FALSE]
#     return(kds)
#   })),nrow=length(mirnas.all),byrow=TRUE)
#   rownames(kds_centered) <- mirnas.all
#   temp_kds <- GetKds("let-7a", "equilibrium", start, stop, sitelist,
#                     combined.input=combined.input, log.residual=log.residual,
#                   nosite=nosite)
#   temp_kd_names <- names(temp_kds)[grep("\\.", names(temp_kds))]
#   colnames(kds_centered) <- c("8mer","6mer", "6mer-m8", temp_kd_names, "None")
#   print(kds_centered)
#   ylims = c(0.5*min(kds_centered), 2*max(kds_centered))
#   plot(1:ncol(kds_centered), kds_centered[1,],
#        type="o",log='y',
#        ylim=ylims,
#        axes=FALSE, 
#        ann=FALSE,
#        pch=19,lwd=2)
#   lines(1:ncol(kds_centered),pch=19,lwd=2,kds_centered[2,],type="o",col="blue")
#   lines(1:ncol(kds_centered),pch=19,lwd=2,kds_centered[3,],type="o",col="red")
#   lines(1:ncol(kds_centered),pch=19,lwd=2,kds_centered[4,],type="o",col="forestgreen")
#   lines(1:ncol(kds_centered),pch=19,lwd=2,kds_centered[5,],type="o",col="purple")
#   print(colnames(kds_centered))
#   axis(1,pos=ylims[1],at=c(1:ncol(kds_centered)),labels=colnames(kds_centered),las=2)
#   axis(2,pos=0.5,at=10^seq(ceiling(log10(min(ylims))-1), floor(log10(max(ylims)))+1))
#   legend("topleft",legend = mirnas.all, col=c("black", "blue", "red", "forestgreen", "purple"), lwd=2, pch=19)

# }



# MakeInputScatter <- function(mirna, experiment, n_constant, ombined.input=TRUE,
#                              log.residual=FALSE, bgoff=FALSE,
#                              vertical.lines=FALSE) {
#   params <- GetSiteKds(mirna, experiment, start, stop, sitelist,
#                    combined.input=combined.input, log.residual=log.residual,
#                    scaled=FALSE)
#   if (combined.input == TRUE) {
#     sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop, sitelist)
#   } else {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
#   }
#   print(sitesXcounts)
#   sitesXcounts <- sitesXcounts[,-1]
#   # Remove any rows with no reads whoatsoever.
#   sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
#   # Remove any rows for which there are no input reads:
#   sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

#   kds <- c(unlist(Logistic(params[1:(nrow(sitesXcounts)-1)], 10)),1)
#   names(kds) <- rownames(sitesXcounts)

#   if (length(site_list) == 0) {
#     site_list <- rownames(sitesXcounts)
#   } else if (class(site_list) == "character") {
#     site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", site_list,".txt"),
#   stringsAsFactors=FALSE)[,1], "None")
#   } else {
#     site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
#   }
#   print(site_list)
#   x <- Norm(sitesXcounts[site_list, 1])
#   y <- Norm(sitesXcounts[site_list, column.plot])
#   xmin <- 0.5 * min(x, y)
#   xmax <- 2 * max(x, y)
#   ymin <- 0.5 * min(x, y)
#   ymax <- 2 * max(x, y)
#   yextension <- (ymax/ymin)
#   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#   xs <- xs[xs >= xmin & xs <= xmax]
#   xmin <- min(xs)
#   xmax <- max(xs)
#   xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
#   ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))
#   ys <- ys[ys >= ymin & ys <= ymax]
#   ymin <- min(ys)
#   ymax <- max(ys)
#   yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))

#   plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.3), ylim=c(ymin, ymax),pch=19, cex = 1.5, ann=FALSE, axes=FALSE,
#      col=kSiteColors[site_list, ])        
#   # Generate tickmarks for axis.
#   axis(1, at=xl,
#        labels=xl*100,
#        pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
#   axis(1, at=xs, labels=FALSE,
#        pos=ymin, lwd=2)
#   # Label the axis at each order of magnitude.

#   axis(2, at=yl,
#        labels=yl*100,
#        pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#   axis(2, at=ys, labels=FALSE,
#        pos=xmin, lwd=2)
#   segments(xmin, xmin, xmax, xmax, lwd = 2, lty = 2)
#   if (vertical.lines == TRUE) {
#     sapply(1:length(x), function(index) {
#       segments(x[index], x[index], x[index], y[index], col=kSiteColors[site_list,][index])
#       })

#   }
#   ago.percent <- as.numeric(colnames(sitesXcounts)[column.plot])/100
#   title(main = paste0(ago.percent*stockago[mirna,"equilibrium"]*100,' pM AGO2-', mirna), font.main=1, cex.main=1.5, line=-2.5, adj=0.1)

#   title(xlab = "Input library (%)", cex.lab=1.5, line=2, adj=0.3)
#   title(ylab = "AGO-bound library (%)", cex.lab=1.5, line=2)
#   print(site_list)
#   print(site_list)
#   legend.names <- rownames(sitesXcounts)[order(kds)]

#   legend.names <- legend.names[which(legend.names %in% site_list)]
#   ordered_list <- legend.names
#   print(ordered_list)
#   centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#                   "11mer-m4.14", "12mer-m4.15", "13mer-m3.15", "11mer-m5.15")
#   centered_rename <- c("Centered-m3.13", "Centered-m3.14", "Centered-m4.14", "Centered.m4-15", "Centered-m3.15", "Centered-m5.15")
#   names(centered_rename) <- centered.sites
#   print(kSiteColors[site_list, ])
#   print(kSiteColors[ordered_list])

#   legend.names <- sapply(legend.names, function(site){
#                         if (site %in% centered.sites) {
#                           return(centered_rename[site])
#                         } else {
#                           return(site)
#                         }
#                        }
#                        )
#   # legend.names <- unique(legend.names)
#   print(legend.names)
#   legend(x=2 * xmax,y=ymax,legend=legend.names, pch=19, col=kSiteColors[ordered_list], cex=1.2, bty="n", ncol = 1)
# }


# MakeWithWithoutProteinScatter <- function(mirna, experiment, start, stop, sitelist) {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
#     colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#       colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#         return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#       }
#     )
#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       k.c.stockago.1 <- 10^params["AGO"]
#       bgs.1 <- 10^params["bg"]

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.2 <- stockago[mirna,experiment]

#       bgs.2 <- 10^params["bg"]

#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
#     print(xmin)
#     print(xmax)
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(c(kds.1, bgs.1, k.c.stockago.1),
#       c(kds.2, bgs.2, k.c.stockago.2),
#       col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
#       pch = 19,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#         segments(xmin, xmin, xmax, xmax, lty = 2)

#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = "Fixed protein parameters (nM)", cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = "Fit protein parameters (nM)", cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=1.2, bty="n", ncol = 1)

    
# }

# MakeSMNBScatter <- function(start, stop, sitelist) {
#     sitesXcounts <- GetSitesXCounts("let-7", "equilibrium", start, stop, sitelist)
#     colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#       colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#         return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#       }
#     )
#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       k.c.stockago.1 <- 10^params["AGO"]
#       bgs.1 <- 10^params["bg"]

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.2 <- stockago[mirna,experiment]

#       bgs.2 <- 10^params["bg"]

#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
#     print(xmin)
#     print(xmax)
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(c(kds.1, bgs.1, k.c.stockago.1),
#       c(kds.2, bgs.2, k.c.stockago.2),
#       col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
#       pch = 19,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#         segments(xmin, xmin, xmax, xmax, lty = 2)

#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = "Fixed protein parameters (nM)", cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = "Fit protein parameters (nM)", cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=1.2, bty="n", ncol = 1)

    
# }

# MakeScatterAcrossWindows <- function(mirna, experiment, lim1, lim2, sitelist, method) {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, lim1, lim1, sitelist)
#     colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#       colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#         return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#       }
#     )
#     if (method == "free_protein") {

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       k.c.stockago.1 <- 10^params["AGO"]
#       bgs.1 <- 10^params["bg"]

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       k.c.stockago.2 <- 10^params["AGO"]
#       bgs.2 <- 10^params["bg"]

#       } else {
#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.1 <- stockago[mirna,experiment]

#       bgs.1 <- 10^params["bg"]


#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.2 <- stockago[mirna,experiment]

#       bgs.2 <- 10^params["bg"]
# }
#     print(kds.1)
#     print(kds.2)
#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
#     print(xmin)
#     print(xmax)
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(c(kds.1, bgs.1, k.c.stockago.1),
#       c(kds.2, bgs.2, k.c.stockago.2),
#       col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
#       pch = 19,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#     segments(xmin, xmin, xmax, xmax, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = paste0(lim1, " parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0(lim2, " parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
# }

# MakeScatterAcrossOrder <- function(mirna, experiment, lim, sitelist1, sitelist2, method) {
#     sitesXcounts <- GetSitesXCounts(mirna, experiment, lim, lim, sitelist1)
#     colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#       colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#         return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#       }
#     )
#     sitesXcounts1 <- sitesXcounts
#         sitesXcounts <- GetSitesXCounts(mirna, experiment, lim, lim, sitelist2)
#     colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#       colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#         return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#       }
#     )
#     sitesXcounts2 <- sitesXcounts
#     print(sitesXcounts2)
#     if (method == "free_protein") {

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
#       k.c.stockago.1 <- 10^params["AGO"]
#       bgs.1 <- 10^params["bg"]

#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg_free_protein.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
#       k.c.stockago.2 <- 10^params["AGO"]
#       bgs.2 <- 10^params["bg"]

#       } else {
#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.1 <- stockago[mirna,experiment]

#       bgs.1 <- 10^params["bg"]


#       params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
#                    "/", experiment, "/kds_with_structure/",
#                    lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg.txt")
#       params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
#       params <- params[nrow(params), - ncol(params)]
#       kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
#       stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                        "/SolveForKds/k_c_stockago.txt"), row.names=1,
#                        header=TRUE, sep="\t")
#       k.c.stockago.2 <- stockago[mirna,experiment]

#       bgs.2 <- 10^params["bg"]
# }
#     names(kds.1) <- rownames(sitesXcounts1)
#     names(kds.2) <- rownames(sitesXcounts2)
#     kds.2 <- kds.2[names(kds.1)]
#     print(kds.1)
#     print(kds.2)
#     xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
#     xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
#     print(xmin)
#     print(xmax)
#     print(bgs.1)
#     print(bgs.2)
#     print(k.c.stockago.1)
#     print(k.c.stockago.2)
#     xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
#     xs <- xs[xs >= xmin & xs <= xmax]
#     xmin <- min(xs)
#     xmax <- max(xs)
#     xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


#     plot(c(kds.1, bgs.1, k.c.stockago.1),
#       c(kds.2, bgs.2, k.c.stockago.2),
#       col = kSiteColors[c(rownames(sitesXcounts1), "bg", "Ago"),],
#       pch = 19,
#       log = 'xy',
#       xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
#       ylim = c(xmin, xmax),
#       axes = FALSE,
#       ann = FALSE)
#     segments(xmin, xmin, xmax, xmax, lty = 2)
#     axis(1, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
#     axis(1, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)
#     # Label the axis at each order of magnitude.

#     axis(2, at=xl,
#          labels=sapply(xl, function(name) {
#            eval(substitute(expression(10^x), list(x=log10(name))))
#          }),
#          pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
#     axis(2, at=xs, labels=FALSE,
#          pos=xmin, lwd=2)

#         title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
#     title(xlab = paste0(sitelist1, " parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
#     title(ylab = paste0(sitelist2, " parameters (nM)"), cex.lab=1.5, line=2)
#     legend.names <- c(rownames(sitesXcounts1), "bg", "Ago")
#     # centered.sites <- c("11mer-m3.13", "12mer-m3.14",
#     #                 "11mer-m4.14", "12mer-m4.15")
#     # centered_index <- sapply(centered.sites, function(site){
#     #                        which(legend.names == site)
#     #                      }
#     #                      )
#     # legend.names[centered_index] <- "Centered"
#     # legend.names <- unique(legend.names)
#     legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
# }






# MakeSiteIterationPlotKineticsSingle <- function(out,method,mirna, num.kds, num.bgs, colors = FALSE) {
#   setEPS()
#   postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/koffs/",
#     mirna,"/", method, ".eps"))
#   par(kPlotParameters)
#   x = seq(1,dim(out)[1],length = max(ncol(out),1000))
#   out <- out[x,]
#   # print(out)
#   probs = out[ ,"-logp"]
#   out.print <- out[ , seq(dim(out)[2] - 1)]
#   out.print.koffs <- Logistic(out[,1:num.kds],max = 10)
#   out.print.setpoint <- 10^out[,num.kds + 1]
#   out.print.bgs <- 10^out.print[,(num.kds + 2) : (num.kds + num.bgs + 1)]

#   out.print <- cbind(out.print.koffs,out.print.setpoint, out.print.bgs)
#   ys <- 10^c(max(floor(log10(min(out.print))), -6), ceiling(log10(max(out.print))))
#   probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])


#   plot(x    = x,
#        y    = 10^probs.scaled,
#        log  = "y",
#        axes = FALSE,
#        ann = FALSE, 
#        type = "l",
#        col = "white",
#        ylim = ys)
#   axis(side = 1,
#        at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
#        pos  = ys[1],
#        tck  = -0.01,
#        labels = FALSE)
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
#   axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
#        labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=1, lwd=2)
#   title(main = mirna, font=1)
#   axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
#        pos=1)
#   if (colors != FALSE) {

#   cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
#   names(cols) <- colnames(out.print)
#   cols["AGO"] <- "grey"

#   sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=cols[name])
#     })
#   } else  {
#     sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=kSiteColors[name,])
#     })
#   }
#   lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
#   dev.off()
# }


# MakeSiteIterationPlotKineticsSingleNoSetPoint <- function(out,method,mirna, num.kds, num.bgs, colors = FALSE) {
#   setEPS()
#   postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/koffs/",
#     mirna,"/", method, ".eps"))
#   par(kPlotParameters)
#   x = seq(1,dim(out)[1],length = max(ncol(out),1000))
#   out <- out[x,]
#   # print(out)
#   probs = out[ ,"-logp"]
#   out.print <- out[ , seq(dim(out)[2] - 1)]
#   out.print.koffs <- Logistic(out[,1:num.kds],max = 10)
#   out.print.bgs <- 10^out.print[,(num.kds + 1) : (num.kds + num.bgs)]

#   out.print <- cbind(out.print.koffs,out.print.bgs)
#   ys <- 10^c(max(floor(log10(min(out.print))), -5), 5)
#   probs.scaled <- probs*0.9 / probs[1] * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])

#   colnames <- c(rep(sapply(colnames(out.print)[1:num.kds], function(name) {
#     trim = unlist(strsplit(name,split = "_"))
#     return(trim[1])
#     }), 3), colnames(out.print)[(3*num.kds + 1):length(colnames(out.print))])
#   names(colnames) <- colnames(out.print)

#   plot(x    = x,
#        y    = 10^probs.scaled,
#        log  = "y",
#        axes = FALSE,
#        ann = FALSE, 
#        type = "l",
#        col = "white",
#        ylim = ys)
#   axis(side = 1,
#        at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
#        pos  = ys[1],
#        tck  = -0.01,
#        labels = FALSE)
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
#   axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
#        labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=1, lwd=2)
#   title(main = mirna, font=1)
#   axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
#        pos=1)
#   if (colors != FALSE) {

#   cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
#   names(cols) <- colnames(out.print)
#   cols["AGO"] <- "grey"

#   sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col="blue")
#     })
#   } else  {
#     sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=kSiteColors[colnames[name],])
#     })
#   }
#   lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
#   dev.off()
#   dev.set(3)
# }

# MakeSiteIterationPlotKineticsDoubleNoSetPoint <- function(out,method,mirna, num.sites, inds_double, num.bgs, colors = FALSE) {
#   setEPS()
#   postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/koffs/",
#     mirna,"/", method, ".eps"))  
#   par(kPlotParameters)
#   x = seq(1,dim(out)[1],length = min(nrow(out),1000))
#   out <- out[x,]
#   probs = out[ ,"-logp"]
#   out.print <- out[ , seq(ncol(out) - 1)]
#   out.print.koffs1 <- 10^out[,1 : num.sites]
#   out.print.koffs2 <- 10^out[,(num.sites + 1) : (2 * num.sites)]
#   out.print.fwd <- 10^out[,(2 * num.sites + 1): (2 * num.sites + length(inds_double))]
#   out.print.rev <- 10^out[,(2 * num.sites + length(inds_double) + 1) : (2 * (num.sites + length(inds_double)))]

#   out.print.bgs <- 10^out.print[,(2 * (num.sites + length(inds_double)) + 1) : (2 * (num.sites + length(inds_double)) + num.bgs)]

#   out.print <- cbind(out.print.koffs1, out.print.koffs2, out.print.fwd, out.print.rev,out.print.bgs)
#   print(dim(out.print))
#   ys <- 10^c(max(floor(log10(min(out.print))), -20), ceiling(log10(max(out.print))))
#   probs.scaled <- probs/(ys[2])
#   colnames <- c(rep(sapply(colnames(out.print)[1:num.sites], function(name) {
#     trim = unlist(strsplit(name,split = "_"))
#     return(trim[1])
#     }), 2),
#     rep(sapply(colnames(out.print)[1:num.sites], function(name) {
#     trim = unlist(strsplit(name,split = "_"))
#     return(trim[1])
#     })[inds.double], 2),
#     colnames(out.print)[(2 * (num.sites + length(inds_double)) + 1):length(colnames(out.print))])
#   print(length(colnames))
#   names(colnames) <- colnames(out.print)
#   out.print <<- out.print
#   plot(x    = x,
#        y    = 10^probs.scaled,
#        log  = "y",
#        axes = FALSE,
#        ann = FALSE, 
#        type = "l",
#        col = "white",
#        ylim = ys)
#   axis(side = 1,
#        at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
#        pos  = ys[1],
#        tck  = -0.01,
#        labels = FALSE)
#   axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
#   axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
#        labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
#          eval(substitute(expression(10^x), list(x=name)))
#        }),
#        pos=1, lwd=2)
#   title(main = mirna, font=1)
#   axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
#        pos=1)
#   print("2821")
#   if (colors != FALSE) {

#   cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
#   names(cols) <- colnames(out.print)
#   cols["AGO"] <- "grey"

#   sapply(colnames(out.print),function(name) {
#     lines(x, out.print[, name], lwd=1, col=cols[name])
#     })
#   } else  {
#     lwds <- c(rep(1, 2*num.sites), rep(1, 2*length(inds.double)), rep(1, num.bgs))
#     ltys <- c(rep(c(1, 2),each=num.sites),rep(c(3, 3),each=length(inds.double)),rep(1, num.bgs))
#     names(lwds) <- colnames(out.print)
#     names(ltys) <- names(lwds)
#     lwds <<- lwds
#     ltys <<- ltys
#     sapply(colnames(out.print),function(name) {

#     lines(x, out.print[, name], lwd=lwds[name],lty = ltys[name] , col=kSiteColors[colnames[name],])
#     })
#   }
#   print("2842")
#   lines(x, probs.scaled, type="l", lwd=2, col="black")
#   dev.off()
#   dev.set(2)
# }


# ModelFunctionDist <- function(pars, c.I.tots, data) {
#   # print(c.totals[1,1])
#   # print(c.I.tots[1])
#   # Split up the parameters into the kd and background parameters.
#   num.kds <- length(c.I.tots)
#   kds  <- 10^pars[1 : num.kds]
#   bgs  <- rep(10^pars[(num.kds + 1)], ncol(data))
#   stock.ago <- 10^pars[num.kds + 2]
#   alpha <- 10^pars[num.kds + 3]
#   beta <- 10^pars[num.kds + 4]
#   k <- 1

#   # Get the bound counts in each exp:
#   c.agos <- sapply(colnames(data), function(x) {
#     as.numeric(x) / 100
#   })
#   c.totals <- matrix(
#                 rep(c.I.tots,ncol(data)), nrow=length(c.I.tots), ncol=ncol(data), byrow=FALSE)


#   n_dist <- 100
#   kds.dist <- c(sapply(1:length(kds), function(kdind) {
#     p.dist <- rbeta(n_dist, shape1=alpha, shape2=beta)
#     k.dist <- (kds[kdind])/(p.dist^k)
#     return(k.dist)
#     }))
#   kds <<- kds
#   kds.dist <<- kds.dist
#   # kds.dist[(length(kds.dist)-n_dist+1):length(kds.dist)] <- kds[length(kds)]
#   c.I.tots.dist <- c(sapply(c.I.tots, function(conc) {
#     rep(conc, n_dist)/n_dist
#   }))
#   c.totals.dist <- matrix(
#                 rep(c.I.tots.dist,ncol(data)), nrow=length(c.I.tots.dist), ncol=ncol(data), byrow=FALSE)
#     c.bounds.dist <- as.matrix(
#     sapply(c.agos, function(percent) {
#       return(GetBoundRNA(kds.dist, c.I.tots.dist, percent * stock.ago))
#     }
#   ))
#   c.bounds.dist.total <- matrix(0, nrow=nrow(data), ncol = ncol(data))
#     for (row in seq(length(c.I.tots))) {
#     rowinds <- (1:n_dist)+(row-1)*n_dist
#     c.bounds.dist.total[row,] <- colSums(c.bounds.dist[rowinds,])
#   }
#   c.frees.dist <- c.totals - c.bounds.dist.total
#   c.bgs.dist <- t(t(c.frees.dist) * bgs / colSums(c.frees.dist))
#   c.all.dist <- c.bounds.dist.total + c.bgs.dist
#   c.final <- data.frame(t(t(c.all.dist) / colSums(c.all.dist) * colSums(data)))

#   return(c.final)
# }

# plotBeta <- function(pars) {
#   num.kds <- length(pars) - 5
#       alpha <- 10^(pars[num.kds + 3])
#     beta  <- 10^(pars[num.kds + 4])
#     x.int <- seq(0, 1, length = 100)

#     plot(x.int, dbeta(x.int, shape1 = alpha, shape2 = beta), type = "l")

# }


# ModelLikelihoodDist <- function(pars, data, input){
#   pars.current <<- pars
#   model <<- ModelFunctionDist(pars, input, data)
#   model.norm <- t(t(model) / colSums(model))
#   data.norm <- t(t(data) / colSums(data))
#   loglikelihood <- -sum(data*log(model.norm))
#   model.R <- model.norm/Norm(input)
#   data.R <- data.norm/Norm(input)
#   loglikelihood <- sum((model.R - data.R)^2)
#   if (tick %% 100 == 0) {
#     par(mfrow=c(1, 2))
#     plot(0, type = "n", log = 'xy', xlim = c(.1,100), ylim=c(0.1, 300))
#     colors.sites <- c(sapply(rownames(data.R)[1:256], GetColorFunction), kSiteColors[rownames(data.R)[257:nrow(data.R)], ])
#     # print(colors.sites)
#     # print(data.R[1:10, ])
#     # print(model.R[1:10, ])
#   sapply(c(seq(1, 256, by =16), 257:nrow(model)), function(row) {
#     points(colnames(data),c(data.R[row,]),col=colors.sites[row], pch = 20)
#     lines(colnames(data),c(model.R[row,]),col=colors.sites[row], pch = 20)
#     num.kds <- nrow(model)
#   })
#     plotBeta(pars)
#     # plot(c(model.R), c(data.R), log = 'xy',col=kSiteColors[rownames(data.R),])
#     print(loglikelihood)
#     print(pars[(length(pars)-3):length(pars)])

#   }

# # }
# tick <<- tick + 1
#   return(loglikelihood)
# }



# CompareDist <- function(pars=GetSiteKds("miR-1", "equilibrium", 5, "paper")) {
#   sitesXcounts <- GetSitesXCounts("miR-1", "equilibrium", 5, "paper")
#   print(sitesXcounts)
#   # Separate site sequences from data file.
#   if (sitelist %in% kKmerSiteLists) {
#     seqs <- rownames(sitesXcounts)
#   } else {
#     seqs <- sitesXcounts[,1]
#     sitesXcounts <- sitesXcounts[,-1]
#   }

#   # Remove any rows with no reads whoatsoever.
#   sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
#   # Remove any rows for which there are no input reads:
#   sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
#   sitesXcounts <<- sitesXcounts
#   sfXc <- GetSiteFlanksXCounts("miR-1", "equilibrium", 5, "paper", "8mer")
#   pars.flanks <- GetFlankKds("miR-1", "equilibrium", 5, "paper", "8mer")
#   sfXc <<- sfXc
#   sitesXcounts <- sitesXcounts[-1,]
#   sitesXcounts <- rbind(sfXc, sitesXcounts)
#   data <- sitesXcounts[,3:7]
#   data0 <- sitesXcounts[,8]
#   names(data0) <- rownames(data)
#   data <- data.matrix(data)
#   data <<- data
#   data0 <<- data0
#   colnames(data) <- colnames(sitesXcounts)[3:7]
#   data.norm <- t(t(data) / colSums(data))
#   I.seq <- sitesXcounts[,2]
#   c.I.tots <- Norm(I.seq)*100
#   c.I.tots <<- c.I.tots
#   c.totals <- matrix(
#                 rep(c.I.tots,ncol(data)), nrow=length(c.I.tots), ncol=ncol(data), byrow=FALSE)
#   colnames(c.totals) <- colnames(data)
#   pars.init <- log10(c(c(pars.flanks$Mean, (pars$Mean)[2:(nrow(pars)-2)])^(0.8), ((pars$Mean)[(nrow(pars)-1)])/10,pars$Mean[nrow(pars)], 1, 2))
  


#   names(pars.init) <- c(rownames(data), "bg", "AGO", "alpha", "beta")
#   print(pars.init)
#   pars.init <<- pars.init
#   opt  <- optim(pars.init, ModelLikelihoodDist, method="BFGS",data=data, input = c.I.tots, control=list(maxit = 100000, parscale = rep(10,length(pars.init))))
#   par.round1 <<- opt$par
#   opt2 <- optim(opt$par, ModelLikelihoodDist, data=data, input = c.I.tots, control=list(maxit = 100000))
#   par.round2 <<- opt2$par
#   opt3 <- optim(opt2$par, ModelLikelihoodDist, data=data, input = c.I.tots, control=list(maxit = 100000))
#   par.round3 <<- opt3$par
#   opt4 <- optim(opt3$par, ModelLikelihoodDist, data=data, input = c.I.tots, control=list(maxit = 100000))
#   par.round4 <<- opt4$par
#   opt5 <- optim(opt4$par, ModelLikelihoodDist, data=data, input = c.I.tots, control=list(maxit = 100000))
#   par.round5 <<- opt5$par

#   print(opt5$par)

#   # c.final <- ModelFunctionDist(c(pars,1), c.I.tots)
#   # data.R <- data.norm/Norm(I.seq)
#   # data.all <- data

#   # print(data.R)

#   # rownames(data)[nrow(data)] <- "None"
#   # print(pars)
#   # kds  <- pars[1 : num.kds]
#   # bgs  <- rep(pars[(num.kds + 1)], ncol(data))
#   # stock.ago <- pars[num.kds + 2]
#   # print(kds)
#   # # Get the bound counts in each exp:
#   # c.agos <- sapply(colnames(data), function(x) {
#   #   as.numeric(x) / 100
#   # })
#   # c.bounds <- as.matrix(
#   #   sapply(c.agos, function(percent) {
#   #     return(GetBoundRNA(kds, c.I.tots, percent * stock.ago))
#   #   }
#   # ))
#   # c.frees <- c.totals - c.bounds
#   # c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
#   # c.all <- c.bounds + c.bgs

#   # # Get the amount of background binding by subtracting the bound from the
#   # # total sites in each exp, normalizing. Must transpose to multiply
#   # # each column.

#   # print(c.final)

#   # e.sd <- 2
#   # n_dist <- 200
#   # kds.dist <- exp(c(sapply(kds, function(kd) {
#   #   e <- log(kd)
#   #   e.dist <- rnorm(n_dist, mean=e, sd=e.sd)
#   #   return(e.dist)
#   #   })))
#   # c.I.tots.dist <- c(sapply(c.I.tots, function(conc) {
#   #   rep(conc, n_dist)/n_dist
#   # }))
#   # c.totals.dist <- matrix(
#   #               rep(c.I.tots.dist,ncol(data)), nrow=length(c.I.tots.dist), ncol=ncol(data), byrow=FALSE)
#   #   c.bounds.dist <- as.matrix(
#   #   sapply(c.agos, function(percent) {
#   #     return(GetBoundRNA(kds.dist, c.I.tots.dist, percent * stock.ago))
#   #   }
#   # ))
#   # c.bounds.dist.total <- c.bounds
#   # for (row in seq(nrow(c.final))) {
#   #   rowinds <- (1:n_dist)+(row-1)*n_dist
#   #   c.bounds.dist.total[row,] <- colSums(c.bounds.dist[rowinds,])
#   # }
#   # c.frees.dist <- c.totals - c.bounds.dist.total
#   # c.bgs.dist <- t(t(c.frees.dist) * bgs / colSums(c.frees.dist))
#   # c.all.dist <- c.bounds.dist.total + c.bgs.dist
#   # c.final.dist <- data.frame(t(t(c.all.dist) / colSums(c.all.dist) * colSums(data)))

#   # c.final.dist.norm <- t(t(c.final.dist) / colSums(c.final.dist))
#   # c.final.R.dist <- c.final.dist.norm/Norm(c.I.tots)
#   #   plot(0, type = "n", log = 'xy', xlim = c(0.003, 1), ylim=c(0.1, 300))
#   # sapply(1:nrow(data.R), function(row) {
#   #   points(c.agos,data.R[row,],col=kSiteColors[rownames(data.R)[row], ], pch = 20)
#   #   lines(c.agos,c.final.R.dist[row,],col=kSiteColors[rownames(data.R)[row], ], pch = 20)

#   # })

# }



# # CompareDist()
# # Iorig.k.m1.6 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kinetics/kmer_counts/I_-3_6.txt")
# # Iorig.k.l7.6 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/I_TGT_-3_6.txt")

# # A4.k.m1.6 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kinetics/kmer_counts/2,1_-3_6.txt")
# # names_temp <- rownames(Iorig.k.m1.6)
# # # I.k.m1.6 <- Norm(A4.k.m1.6[,1,drop=FALSE])/Norm(Iorig.k.m1.6[,1,drop=FALSE])
# # # I.k.l7.6 <- Norm(A4.k.m1.6[,1,drop=FALSE])/Norm(Iorig.k.l7.6[,1,drop=FALSE])

# # I.k.m1.6 <- Norm(Iorig.k.l7.6[,1,drop=FALSE])*10000
# # I.k.l7.6 <- Norm(A4.k.m1.6[,1,drop=FALSE])*10000


# # rownames(I.k.m1.6) <- names_temp
# # rownames(I.k.l7.6) <- rownames(I.k.m1.6)
# # spike124 <- "TCTCATCTACCTCCCGGTTTTAATGAATAGTGCCTTAGAC"
# # spike430 <- "TCTCATCTACCTCCCGGTTTTAATGAATAAGCACTTAGAC"
# # spike427 <- "TCTCATCTACCTCCCGGTTTTAATGAATAGCACTTTCGAC"
# # spike14 <- "TCTCATCTACCTCCCGGTTTTAATGAATATCGCTCCCGAC"
# # spike155 <- "TCTCATCTACCTCCCGGTTTTAATGAATAAGCATTAAGAC"

# # spike124rc <- "GTCTAAGGCACTATTCATTAAAACCGGGAGGTAGATGAGA"
# # spike430rc <- "GTCTAAGTGCTTATTCATTAAAACCGGGAGGTAGATGAGA"
# # spike427rc <- "GTCGAAAGTGCTATTCATTAAAACCGGGAGGTAGATGAGA"
# # spike14rc  <- "GTCGGGAGCGATATTCATTAAAACCGGGAGGTAGATGAGA"
# # spike155rc <- "GTCTTAATGCTTATTCATTAAAACCGGGAGGTAGATGAGA"



# # constant_5pw1 <- "GATCGTCGGACTGTAGAACTCCGCCCC"
# # constant_5pw2 <- "GATCGTCGGACTGTAGAACCCTGCCCC"
# # constant_5pw3 <- "GATCGTCGGACTGTAGAACCCCGCCCC"

# # constant_5p   <- "GATCGTCGGACTGTAGAACTCTGCCCC"

# # constant_3p <- "CAAGCAGAAGACGGCATACGA"

# # inds_5pw1 <- sapply(rownames(I.k.m1.6), grepl, x=constant_5pw1,fixed=TRUE)
# # inds_5pw2 <- sapply(rownames(I.k.m1.6), grepl, x=constant_5pw2,fixed=TRUE)
# # inds_5pw3 <- sapply(rownames(I.k.m1.6), grepl, x=constant_5pw3,fixed=TRUE)

# # inds_5p <- sapply(rownames(I.k.m1.6), grepl, x=constant_5p,fixed=TRUE)
# # inds_3p <- sapply(rownames(I.k.m1.6), grepl, x=constant_3p, fixed=TRUE)

# # inds_s124   <- sapply(rownames(I.k.m1.6), grepl, x=spike124, fixed=TRUE)
# # inds_s124rc <- sapply(rownames(I.k.m1.6), grepl, x=spike124rc, fixed=TRUE)

# # inds_s430     <- sapply(rownames(I.k.m1.6), grepl, x=spike430, fixed=TRUE)
# # inds_s430rc   <- sapply(rownames(I.k.m1.6), grepl, x=spike430rc, fixed=TRUE)

# # inds_s427     <- sapply(rownames(I.k.m1.6), grepl, x=spike427, fixed=TRUE)
# # inds_s427rc   <- sapply(rownames(I.k.m1.6), grepl, x=spike427rc, fixed=TRUE)

# # inds_s14      <- sapply(rownames(I.k.m1.6), grepl, x=spike14, fixed=TRUE)
# # inds_s14rc    <- sapply(rownames(I.k.m1.6), grepl, x=spike14rc, fixed=TRUE)

# # inds_s155     <- sapply(rownames(I.k.m1.6), grepl, x=spike155, fixed=TRUE)
# # inds_s155rc   <- sapply(rownames(I.k.m1.6), grepl, x=spike155rc, fixed=TRUE)



# # cols <- rep("gray", length=nrow(I.k.m1.6))
# # cols[inds_5pw1] <- "cyan"
# # cols[inds_5pw2] <- "cyan"
# # cols[inds_5pw3] <- "cyan"

# # cols[inds_5p] <- "blue"
# # cols[inds_3p] <- "red"
# # plot(I.k.m1.6[,1], I.k.l7.6[,1],log='xy',col=cols,xlim = c(0.1, 500), ylim = c(0.1, 500))
# # x_range <- 10^(seq(-2, 5))
# # lines(x_range,x_range,lty=2)
# # lines(x_range, 2*x_range,lty=2, col="gray")
# # lines(x_range, 0.5*x_range, lty = 2, col = "gray")
# # points(I.k.m1.6[inds_5p,1], I.k.l7.6[inds_5p,1],log='xy',col="blue")
# # points(I.k.m1.6[inds_3p,1], I.k.l7.6[inds_3p,1],log='xy',col="red")
# # points(I.k.m1.6[inds_s430,1], I.k.l7.6[inds_s430,1],log='xy',col="forestgreen")
# # points(I.k.m1.6[inds_s430rc,1], I.k.l7.6[inds_s430rc,1],log='xy',col="forestgreen")
# # points(I.k.m1.6[inds_s427,1], I.k.l7.6[inds_s427,1],log='xy',col="green")
# # points(I.k.m1.6[inds_s427rc,1], I.k.l7.6[inds_s427rc,1],log='xy',col="green")
# # points(I.k.m1.6[inds_s14,1], I.k.l7.6[inds_s14,1],log='xy',col="limegreen")
# # points(I.k.m1.6[inds_s14rc,1], I.k.l7.6[inds_s14rc,1],log='xy',col="limegreen")
# # points(I.k.m1.6[inds_s124,1], I.k.l7.6[inds_s124,1],log='xy',col="purple")
# # points(I.k.m1.6[inds_s124rc,1], I.k.l7.6[inds_s124rc,1],log='xy',col="purple")
# # points(I.k.m1.6[inds_s155,1], I.k.l7.6[inds_s155,1],log='xy',col="violet")
# # points(I.k.m1.6[inds_s155rc,1], I.k.l7.6[inds_s155rc,1],log='xy',col="violet")


# # for (ind_row in seq(nrow(I.k.m1.6[inds_5p,]))) {
# #   pos = regexpr(rownames(I.k.m1.6[inds_5p,])[ind_row],constant_5p)
# #   print("pos")
# #   print(pos)[1]
# #   print("row:")
# #   print(I.k.m1.6[inds_5p,][ind_row,])
# #   text(I.k.m1.6[inds_5p,][ind_row,1]*rnorm(1,1,0.05),I.k.l7.6[inds_5p,][ind_row,1]+rnorm(1,1,0.05),pos)
# # }
# # break
# # identify(I.k.m1.6[,1], I.k.l7.6[,1],labels=rownames(I.k.m1.6))
# # # I_TGT.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/I_TGT_-3_6.txt")
# # # I_ACA.kpp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/I_ACA_-3_6.txt")

# # # t0.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/0_-3_6.txt")
# # # t3.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/3_-3_6.txt")
# # # t5.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/5_-3_6.txt")
# # # t10.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/10_-3_6.txt")
# # # t30.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/30_-3_6.txt")
# # # t90.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/90_-3_6.txt")
# # # t300.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/300_-3_6.txt")
# # # t900.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/900_-3_6.txt")
# # # t2700.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/2700_-3_6.txt")
# # # e.kp.m1 <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/kin_pilot/kmer_counts/Equil_-3_6.txt")

# # break

