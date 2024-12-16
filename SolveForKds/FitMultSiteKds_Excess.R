dinucs <- read.table(file=sprintf("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/%s/dinucleotide_counts/I_5.txt",
                                  mirna, experiment)

rownames(dinucs) <- sort(kDinucs)

GetSequenceKmers <- function(x, kmer_l) {
  s1 <- substr(x, 1, nchar(x) - 1)
  s2 <- substr(x, 2, nchar(x))
  sapply(seq(1, nchar(x) - kmer_l + 1), function(l) {
    substr(x, l, l + kmer_l - 1)
  })
}

GetKmerIndex <- function(x) {
  nt_inds <- c("A"=0, "C"=1, "G"=2, "T"=3)
  bases <- 4^(seq(nchar(x)-1, 0))
  values <- nt_inds[unlist(strsplit(x, split=""))]
  return(sum(values*bases) + 1)
}

MakeKmerProbMatrix <- function(mirna, experiment, condition, n_constant,
                               kmer_length, constant=TRUE) {
  kmers_raw <- SubfunctionCall(GetPositionalKmers)
  kmers_init <- MatNorm(kmers_raw)
  kmers_conditional <- sapply(seq(0, nrow(kmers_raw)/4 - 1), function(pre_ind) {
    mat_inds <- 4*pre_ind + c(1, 2, 3, 4)
    totals <- colSums(kmers_raw[mat_inds, ])
    totals_mat <- matrix(totals, nrow=4, ncol=length(totals), byrow=TRUE)
    rownames(totals_mat) <- rownames(kmers_raw)[mat_inds]
    colnames(totals_mat) <- colnames(kmers_raw)
    list(totals_mat)
    })
  kmers_conditional <- do.call(rbind, kmers_conditional)
  inds_nonzero <- which(kmers_conditional != 0)
  kmers_conditional[inds_nonzero] <- (unlist(kmers_raw)[inds_nonzero])/(c(kmers_conditional)[inds_nonzero])
  if (!(constant)) {
    kmers_init <- kmers_init[, 6:(6 + 37 - 1)]
    kmers_conditional <- kmers_conditional[, 6:(6 + 37 - 1)]

  }
  list(kmers_init, kmers_conditional)
}

NextStartPosition <- function(mirna, site1, site2) {
  seq1_start <- GetSiteSeq(mirna, site1)
  seq2_start <- GetSiteSeq(mirna, site2)
  seq1 <- substr(seq1_start, 2, nchar(seq1_start))
  seq2 <- substr(seq2_start, 1, nchar(seq2_start) - 1)
  shift <- 1
  while(seq1 != seq2) {
    seq1 <- substr(seq1, 2, nchar(seq1))
    seq2 <- substr(seq2, 1, nchar(seq2) - 1)
    shift <- shift + 1
  }
  space_str = paste0(rep(" ", length=shift), collapse="")
  return(shift) 
}

GetSequenceOverlap <- function(seq1_start, seq2_start) {
  seq1 <- substr(seq1_start, 2, nchar(seq1_start))
  seq2 <- substr(seq2_start, 1, nchar(seq2_start) - 1)
  end <- substr(seq2_start, nchar(seq2_start), nchar(seq2_start))
  while(seq1 != seq2 & (seq1 != "" | seq2 != "")) {
    end <- paste0(substr(seq2, nchar(seq2), nchar(seq2)), end)
    seq1 <- substr(seq1, 2, nchar(seq1))
    seq2 <- substr(seq2, 1, nchar(seq2) - 1)
  }
  overlap <- paste0(seq1_start, end)
  if (overlap == paste0(seq1_start, seq2_start)) {
    return("")
  } else {
    return(overlap)
  }
}

GetKmerProbability <- function(sequence, kmer_l, pos=1, prior=NULL,
                               constant=TRUE) {
  kmer_prob <- MakeKmerProbMatrix(mirna, experiment, "I", n_constant, kmer_l,
                                  constant=constant)
  kmer_prob <<- kmer_prob
  seq_kmers <- GetSequenceKmers(sequence, kmer_l)
  if (nchar(sequence) < kmer_l) {
    # Grep to get all the matching kmers within the n-nucleotide frequencies
    prob <- sum(kmer_prob[[1]][grep(sprintf("^%s", sequence),
                                    rownames(kmer_prob[[1]]), perl=TRUE), pos])
  } else {
    # Make the dinucleotide frequency tables:
    seq_inds <- sapply(seq_kmers, GetKmerIndex)
    # The the frequency of the first dinucleotide:
    if (length(prior) != 0) {
      prior_nuc = prior[1]
      prior_pos = as.numeric(prior[2])
      N_pos <- max(0, pos - prior_pos - 1)
      overlap <- -1*min(0, pos - prior_pos - 1)
      if (overlap) {
        start_prob <- kmer_prob[[2]][seq_inds[1], pos]    
      } else if (N_pos == 0) {
        prior_seq = paste0(prior_nuc, substr(sequence, 1, 1))
        prior_seq_ind <- GetKmerIndex(prior_seq)
        start_prob <- (kmer_prob[[2]][prior_seq, prior_pos] *
                       kmer_prob[[2]][seq_inds[1], pos])
      } else {
        left_str <- grep(sprintf("^%s", prior_nuc), rownames(kmer_prob[[2]]),
                         perl=TRUE)
        left <- kmer_prob[[2]][left_str, prior_pos]
        right_str <- grep(sprintf("%s$", substr(sequence, 1, 1)),
                          rownames(kmer_prob[[2]]), perl=TRUE)
        right <- kmer_prob[[2]][right_str, pos - 1]
        if (N_pos > 1) {
          middle_mat <- kmer_prob[[2]][,(prior_pos + 1):(pos - 2), drop=FALSE]
          for (col in seq(ncol(middle_mat))) {
            x <- middle_mat[, col]
            dim(x) <- c(4^(kmer_l - 1), 4^(kmer_l - 1))
            x <- t(x)
            left <- left %*% x
          }
          start_prob <- (left %*% right) * kmer_prob[[2]][seq_inds[1], pos]
        } else { # N = 1
          start_prob <- (left %*% right) * kmer_prob[[2]][seq_inds[1], pos]
        }
      }
    } else {
      start_prob <- kmer_prob[[1]][seq_inds[1], pos]
    }
    # The the frequency of the rest:
    if (nchar(sequence) == kmer_l) {
      further_prob <- 1
    } else {
      further_prob <- sapply(2:length(seq_inds), function(ind) {
                             kmer_prob[[2]][seq_inds[ind], pos + ind - 1]
                           })
    }
  prob <- exp(log(start_prob) + sum(log(further_prob)))
}
  prob
}


ProbAtLeast1 <- function(site, kmer_l) {
  sequence <- GetSiteSeq(mirna, site)
  total_positions <- 37 + n_constant*2 - nchar(sequence) + 1
  all_probs <- sapply(1:total_positions, GetKmerProbability, sequence=sequence,
                      kmer_l=kmer_l)
  return(1-exp(sum(log(1 - all_probs))))
}

ProbExactly1 <- function(site, kmer_l, raw=FALSE) {
  if (!raw) {
    sequence <- GetSiteSeq(mirna, site)  
  } else {
    sequence <- site
  }
  total_positions <- 37 + n_constant*2 - nchar(sequence) + 1
  spacing <- NextStartPosition(mirna, site, site)
  all_probs <- sapply(1:total_positions, function(pos) {
    occluded_spaces <- seq(pos - spacing + 1, pos + spacing - 1)
    occluded_spaces <- occluded_spaces[occluded_spaces > 0]
    occluded_spaces <- occluded_spaces[occluded_spaces <= max(total_positions)]
    prob_pos <- GetKmerProbability(sequence, kmer_l, pos=pos, constant=TRUE)
    prior <- c(substr(sequence, nchar(sequence), nchar(sequence)),
               nchar(sequence))

    # random_string = paste0(rep("n", 37), collapse="")
    # constant_string = paste0(rep("-", 5), collapse="")

    # string_space_5p = paste0(rep("_", pos - 1), collapse="")
    # string_space_3p = paste0(rep("_", 37 + n_constant*2 - nchar(sequence) - pos + 1), collapse="")
    # message(paste0(5, 4, 3, 2, 1, random_string, 1, 2, 3, 4, 5))

    # message(sprintf("%s%s%s %s", string_space_5p, sequence, string_space_3p, prob_pos))
    remaining_positions <- setdiff(seq(total_positions), occluded_spaces)
    non_probs <- sapply(remaining_positions, function(pos1) {
                        1 - SubfunctionCall(GetKmerProbability, prior=prior, pos=pos1, constant=TRUE)})
    for (i in seq(length(non_probs))) {
      remain_pos <- remaining_positions[i]
      non_prob <- non_probs[i]
      # print(remain_pos)
      # print("post repeat unit:")
      # print(total_positions - spacing - remain_pos + 1)
      prespace = sprintf(paste0(rep(".", remain_pos - 1), collapse=""))
      postspace = sprintf(paste0(rep(".",
                                     total_positions - remain_pos), collapse=""))
      # message(paste0(prespace, tolower(sequence), postspace, " ", sprintf("%s", non_prob)))
      # message(paste0(paste0(rep(" ", i), collapse=""), sprintf("%s", non_prob)))
    }
    # message(sprintf("%s%s%s %s", string_space_5p, sequence, string_space_3p, prob_pos))

    final_prob <- exp(sum(log(c(prob_pos, non_prob))))
    # print(final_prob)
    final_prob
  })
  sum(all_probs)
}

ProbExactly1_Simple <- function(site, kmer_l, raw=FALSE) {
  if (!raw) {
    sequence <- GetSiteSeq(mirna, site)  
  } else {
    sequence <- site
  }
  total_positions <- 37 + n_constant*2 - nchar(sequence) + 1
  sum(sapply(1:total_positions, function(pos) {
    SubfunctionCall(GetKmerProbability, prior=NULL, pos=pos, constant=TRUE)
  }))
}

PlotSingleSiteFrequencies <- function() {
  graphics.off()
  dev.new(xpos=20, ypos=20, height=5, width=5)
  par(kPlotParameters)
  xmin <- 0
  xmax <- ncol(single_fracs)
  ymin <- 0
  ymax <- 1
  BlankPlot()
  x <- 1:ncol(single_fracs)
  print(x)
  textpos <- seq(1, 7, length.out=nrow(single_fracs))
  for (i in seq(nrow(single_fracs))) {
    lines(x, single_fracs[i, ], col=kSiteColors[rownames(single_fracs)[i]])
    text(textpos[i], single_fracs[i,5], labels=rownames(single_fracs)[i], col=kSiteColors[rownames(single_fracs)[i]])
  }
  xmin <- 1
  AddLinearAxis(1, tick.space=1, label.space=1, label="Sample", alt_lab=colnames(single_fracs))
  AddLinearAxis(2, tick.space=0.2, label.space=0.2, label="Fraction of reads singly bound")
  xy <- GetPlotFractionalCoords(fx=0.9, fy=0.05)
  text(xy[1], xy[2], labels=mirna, adj=1)
  # dev.copy2pdf(file=sprintf("figures/%s_fraction_of_single_sites.pdf", mirna))
}
