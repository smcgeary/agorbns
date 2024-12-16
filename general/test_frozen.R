


# mirna1 <- "let-7a"
# mirna2 <- "let-7a"

experiment1 <- "equilibrium_mmseed_nb"
experiment2 <- experiment1

Determine8merReferenceSiteType <- function(mirna, kmer_vec) {
  # Make the vector of letters in the 8mer site type.
  site_8mer <- StrRev(GetSiteSeq(mirna, "8mer"))
  site_8mer_vec <- unlist(strsplit(site_8mer, split=""))
  # Apply a loop over the kmers.
  site_vec <- sapply(kmer_vec, function(kmer) {
    # Condition in case the site is "None"
    if (kmer == "None") {
      return(kmer)
    # Condition in which the kmer is actually a kmer.
    } else {
      suffix <- ""
      kmer_vec <- unlist(strsplit(StrRev(kmer), split=""))
      for (i in 1:8) {
        if (site_8mer_vec[i] != kmer_vec[i]) {
          suffix <- sprintf("%smm%s%s", suffix, kmer_vec[i], i)
        } 
      }
      if (suffix != "") {
        suffix <- sprintf("-%s", suffix)
      }
      out <- sprintf("8mer%s", suffix)
      return(out)
    }
  })
  return(site_vec)
}

counts1 <- SubfunctionCall(GetProgrammedKmerFrequencies, mirna=mirna1,
                           experiment=experiment1, condition="I")
counts2 <- SubfunctionCall(GetProgrammedKmerFrequencies, mirna=mirna2,
                           experiment=experiment2, condition="I")
# Determine the names of each of the possible 8mer sites (1), 7mer-m8 sites 
# (3), 7mer-A1 sites (3), 6mer sites (9), and 8mer-mm sites (18), for both
# miRNAs.



sites_alt <- Determine8merReferenceSiteType(mirna2, rownames(counts1))
sites_alt <- sites_alt[-length(sites_alt)]

# GetOverlapMatrix <- function(sites_alt) {
#   sites_single_mismatch <- sites_alt[5:(18 + 4)]
#   # print(sites_single_mismatch)
#   mat_nonzero <- sapply(sites_single_mismatch, function(site) {
#     # print(site)
#     suffix <- unlist(strsplit(site, split="-"))[2]
#     indeces <- grep(suffix, sites_alt)
#     # print(indeces)
#     indeces[1] <- 1
#     row <- rep(0, length(sites_alt))
#     row[indeces] <- 1
#     row 
#   })
#   full_matrix <- matrix(rep(0, length(sites_alt)^2), nrow=length(sites_alt),
#                         ncol=length(sites_alt))
#   full_matrix[, 5:(18 + 4)] <- mat_nonzero
#   return(full_matrix)
# }

# GetOverlapMatrix2 <- function(mirna, sites_alt, nucs_dict) {
#   # First make the final matrix, pre-allocated with all zeros.
#   full_matrix <- matrix(rep(0, length(sites_alt)^2), nrow=length(sites_alt),
#                       ncol=length(sites_alt))
#   temp_col <- rep(0, length(sites_alt))
#   # Make the list of miRNA letters.
#   site_8mer_vec <- unlist(strsplit(StrRev(GetSiteSeq(mirna, "8mer")), split=""))
#   sites_single_mismatch <- sites_alt[2:25]
#   col_8mer <- temp_col
#   sites_8mer1mm_pos <- gsub("8mer-mm", replacement="", x=sites_single_mismatch)
#   sites_8mer1mm_pos <- as.integer(substring(sites_8mer1mm_pos, first=2))
#   col_8mer[2:25] <- nucs_dict[site_8mer_vec[sites_8mer1mm_pos]]
#   # print(col_8mer)
#   full_matrix[, 1] <- col_8mer
#   mat_nonzero <- sapply(sites_single_mismatch, function(site) {
#     suffix <- unlist(strsplit(site, split="-"))[2]
#     mm <- substring(suffix, first=3, last=3)
#     pos <- as.integer(substring(suffix, first=4))
#     indeces <- grep(suffix, sites_alt)
#     # Series of three string operations intended to determine the identify of
#     # the mutated nucleotide.
#     values <- grep(suffix, sites_alt, value=TRUE)
#     values <- gsub("8mer-", replacement="", x=values)[-1]
#     values <- substring(
#       gsub(suffix, replacement="", x=values), first=3, last=3
#     )
#     values <- c(site_8mer_vec[pos], values)
#     # This needs to be pre-appended to the identity of the nucleotide that the
#     # mismatch itself is, since mutation back makes it into the 8mer.
#     indeces[1] <- 1
#     # row <- rep(0, length(sites_alt))
#     temp_col[indeces] <- nucs_dict[values]
#     temp_col
#   })

#   full_matrix[, 2:25] <- mat_nonzero

#   double_matrix <- sapply(26:length(sites_alt), function(i_col) {
#     site <- sites_alt[i_col]
#     site <- gsub("8mer-", replacement="", x=site)
#     double_site <- unlist(strsplit(site, split="mm"))[-1]
#     pos <- as.integer(rev(substring(double_site, first=2)))
#     vals <- nucs_dict[site_8mer_vec[pos]]
#     inds <- sapply(double_site, function(site) {
#       which(sites_alt == sprintf("8mer-mm%s", site))
#     })
#     temp_col[inds] <- vals
#     temp_col
#   })
#   full_matrix[, 26:length(sites_alt)] <- double_matrix
#   return(full_matrix)
# }

GetSiteMutationMatrix <- function(mirna) {
  # First make the final matrix, pre-allocated with all zeros.
  counts <- GetProgrammedKmerFrequencies(mirna, "equilibrium_mmseed_nb",
                                         condition="I")

  counts <- counts[-nrow(counts), , drop=FALSE]

  sites <- Determine8merReferenceSiteType(mirna, rownames(counts)[])
  full_matrix <- matrix(rep("None", (length(sites) + 1)^2),
                        nrow=length(sites) + 1,
                        ncol=length(sites) + 1)
  temp_col <- rep("None", length(sites) + 1)
  # Make the list of miRNA letters.
  site_8mer_vec <- unlist(strsplit(StrRev(GetSiteSeq(mirna, "8mer")), split=""))
  # Make the first column, which describes mutation from the 8mer site
  # to one of the 18 1-nt mismatches.
  sites_single_mismatch <- sites[2:25]
  col_8mer <- temp_col
  s8mer1mm <- gsub("8mer-mm", replacement="", x=sites_single_mismatch)
  pos <- as.integer(substring(s8mer1mm, first=2))
  nuc_start <- site_8mer_vec[pos]
  nuc_end <- substring(s8mer1mm, first=1, last=1)
  nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
  col_8mer[2:25] <- nuc_combined
  full_matrix[, 1] <- col_8mer
  # Make the next 24 column, which describes mutation from one of the 8mer-1mm
  # sites into either the 8mer (the first entry) or 
  # to one of the 18 1-nt mismatches.
  cols_1mm <- sapply(sites_single_mismatch, function(site) {
    # Get the 1mm suffix, allocate the mismatch nucleotide and pos.
    suffix <- unlist(strsplit(site, split="-"))[2]
    mm1_nuc <- substring(suffix, first=3, last=3)
    mm1_pos <- as.integer(substring(suffix, first=4))
    # Get the indeces of which column entries will used. Note: The first is
    # incorrect and will be changed to `1`, because each 1mm-site can mutate to
    # the 8mer.
    inds <- grep(suffix, sites)
    inds[1] <- 1
    # Series of three string operations intended to determine the identify of
    # the mutated nucleotide.
    mm2_strings <- gsub("8mer-", replacement="",
                        x=grep(suffix, sites, value=TRUE))[-1]
    # Remove the starting mismatch from these strings to yield the second
    # mismatch
    mm2 <- gsub(suffix, replacement="", x=mm2_strings)
    mm2_pos <- as.integer(substring(mm2, first=4, last=4))
    mm2_nuc <- substring(mm2, first=3, last=3)

    nuc_start <- c(mm1_nuc, site_8mer_vec[mm2_pos])
    nuc_end <- c(site_8mer_vec[mm1_pos], mm2_nuc)
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    temp_col[inds] <- nuc_combined
    # Get the mutations to 1mm sites at the same position (2 possible).
    # gsub_string <- sprintf("^8mer-mm.%s$", mm1_pos)
    # # The starting nucleotide is the mismatch nucleotide.
    # nuc_start <- rep(mm1_nuc, 2)
    # # Get the index of the other two mutants at this position.
    # inds <- setdiff(grep(gsub_string, sites, perl=TRUE),
    #                 c(grep(sprintf("^%s$", site), sites, perl=TRUE)))
    # nuc_end <- substring(sites[inds], first=8, last=8)
    # nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    # temp_col[inds] <- nuc_combined
    temp_col
  })
  # Allocate these columns
  full_matrix[, 2:25] <- cols_1mm
  # Now determine the indeces for the remaining double sites mutating to single
  # sites.
  double_matrix <- sapply(26:length(sites), function(i_col) {
    # First deal with reversion
    mm2 <- gsub("8mer-", replacement="", x=sites[i_col])
    mm2_vec <- unlist(strsplit(mm2, split="mm"))[-1]
    mm2_pos <- as.integer(rev(substring(mm2_vec, first=2)))
    nuc_start <- rev(substring(mm2_vec, first=1, last=1))
    nuc_end <- site_8mer_vec[mm2_pos]
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    inds <- sapply(mm2_vec, function(site) {
      which(sites == sprintf("8mer-mm%s", site))
    })
    temp_col[inds] <- nuc_combined
    # Now deal with mutation to one of the other two for each.
    gsub_1 <- sprintf("8mer-mm.%smm%s", mm2_pos[2], mm2_vec[2])
    gsub_2 <- sprintf("8mer-mm%smm.%s", mm2_vec[1], mm2_pos[1])
    nuc_start <- rep(c(substring(mm2, first=3, last=3),
                       substring(mm2, first=7, last=7)), each=2)
    inds <- setdiff(c(grep(gsub_1, sites),
                      grep(gsub_2, sites)), c(i_col))
    other_mm_sites <- gsub("8mer-", replacement="", x=sites[inds])
    nuc_end <- c(substring(other_mm_sites[1:2], first=3, last=3),
                 substring(other_mm_sites[3:4], first=7, last=7))
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    temp_col[inds] <- nuc_combined
    temp_col
  })
  full_matrix[, 26:length(sites)] <- double_matrix
  print(length(sites))
  print(dim(full_matrix))
  rownames(full_matrix) <- c(sites, "triple_mm")
  colnames(full_matrix) <- c(sites, "triple_mm")
  return(full_matrix)
}


M_l7 <- GetSiteMutationMatrix("let-7a")
M_m1 <- GetSiteMutationMatrix("miR-1")
M_m155 <- GetSiteMutationMatrix("miR-155")
M_list <- list(`let-7a`=M_l7, `miR-1`=M_m1, `miR-155`=M_m155)

print(M_l7[1:25, 1:25])

pc <- 1e-8

# M2 <- GetOverlapMatrix2(mirna2, sites_alt, exp(pars_init))

pars_dinuc_init <- rep(-5, 16)
names(pars_dinuc_init) <- kDinucs


pars_dinuc_init <- pars_dinuc_init[c(-1, -4, -13, -16)]

print(pars_dinuc_init)

tick <- 0

# M3 <- GetOverlapMatrix3(mirna2, sites_alt)
M3 <- M_m155
FitLibraries <- function(pars, mirna=mirna2) {
  pars_extended <- c(exp(pars), "None"=0)
  M_use <- M_list[[mirna2]]
  M <- matrix(pars_extended[M_use], nrow=nrow(M_use))
  M[nrow(M3), 26:(nrow(M3) - 1)] <- 18*mean(exp(pars))
  rownames(M) <- rownames(M3)
  colnames(M) <- colnames(M3)
  start_counts <- rep(0, nrow(counts1))
  names(start_counts) <- rownames(counts1)[length(start_counts)]
  counts <- counts2[1:length(start_counts), 1, drop=FALSE]
  if (mirna == "miR-155") {
    inds <- grep("^8mer-mm.1mm.8$", rownames(M), perl=TRUE)
    inds_use <- c(5:22, inds)
  } else {
    inds_use <- 5:22
  }
  start_counts[inds_use] <- counts[inds_use, ]
  model <- (diag(1 - colSums(M)) + M)%*%(diag(1 - colSums(M)) + M)%*%start_counts + pc
  counts <- counts + pc
  model <- model/sum(model)
  cost <- sum(-1*counts[-nrow(M), 1]*log(model[-nrow(M)]))
  tick <<- tick + 1
  if (tick%%100 == 0) {
    plot(model[-nrow(M)], counts[-nrow(M), 1], log='xy',
       xlim=c(1e-8, 1), ylim=c(1e-8, 1))
    segments(1e-8, 1e-8, x1=1, y1=1, lty=2)
    print(cost)
  }
  cost
}

opt <- optim(par=pars_dinuc_init, fn=FitLibraries, method="BFGS")
opt2 <- optim(par=opt$par, fn=FitLibraries, method="BFGS")
opt3 <- optim(par=opt2$par, fn=FitLibraries, method="BFGS")
opt4 <- optim(par=opt3$par, fn=FitLibraries, method="BFGS")