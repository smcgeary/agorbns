


GetThreePrimeKdMatrix <- function(mirna, register, kmer_len) {
	# This gets Kd matrices for full complementary kmers, as calculated by
	# Namita.
	full_path <- sprintf("/%s/%s/matdfchange%smer%s_reg%s.txt", base_path,
	                     mirna_path[mirna], kmer_len, mirna_path[mirna], register)
	print(full_path)
	dG_table <- read.table(full_path, sep=",", row.names=1, header=TRUE)
	colnames(dG_table) <- gsub("^(X)", colnames(dG_table), replace="")
	dG_table
}


kMirnasThrP <- c("let-7a", "miR-1", "miR-155")
kLibsThrP <- c("mix2_rep", "mix2", "mix2_3_3")

kKDTypesThrP <- c("kds_site_counts_streamline3p_longloop_pyr-pyrsWCs",
                  "kds_site_counts_streamline3p_longloop_pur-pursWCs",
                  "kds_site_counts_streamline3p_longloop_ACsWCs",
                  "kds_site_counts_streamline3p_longloop_GUsWCs",
                  "kds_site_counts_streamline3p_longloop_BsWCs")


# Dictionary for converting the base names in Namita's Kd tables.
kSiteNameConversion <- list(
  `let-7a`=c(
		`8mer-1.8`="8mer", `7mer-2.8`="7mer_m8", `7mer-1.7`="7mer_A1",
		`6mer-2.7`="6mer", `6mer-3.8`="6mer_m8", `6mer-1.6`="6mer_A1",
		`8mer-1.8_mm.AA.7`="8mer_mmA7", `8mer-1.8_mm.AC.7`="8mer_mmC7",
		`8mer-1.8_mm.TT.6`="8mer_mmT6", `8mer-1.8_mm.TC.6`="8mer_mmC6",
		`8mer-1.8_mm.GT.2`="8mer_w2", `8mer-1.8_mm.AC.3`="8mer_mmC3",
		`8mer-1.8_mm.AG.7`="8mer_mmG7", `8mer-1.8_mm.TG.6`="8mer_w6",
		`8mer-1.8_mm.AA.3`="8mer_mmA3", `8mer-1.8_mm.GA.2`="8mer_mmA2",
		`8mer-1.8_mm.GG.2`="8mer_mmG2", `8mer-1.8_mm.GT.5`="8mer_w5",
		`8mer-1.8_mm.GT.4`="8mer_w4", `8mer-1.8_mm.GA.4`="8mer_mmA4",
		`8mer-1.8_mm.GA.5`="8mer_mmA5", `8mer-1.8_mm.AG.3`="8mer_mmG3",
		`8mer-1.8_mm.GG.5`="8mer_mmG5", `8mer-1.8_mm.GG.4`="8mer_mmG4",
		`None`="None", `bg`="bg", `AGO`="AGO")
)




GetKdPath <- function(mirna, experiment, library_type_ind, kd_type_ind) {
	# This gets Kd matrices for full complementary kmers, as calculated by
	# Namita.
	library_type <- kLibsThrP[library_type_ind]
	kd_type <- kKDTypesThrP[kd_type_ind]
	directory <- sprintf("/lab/solexa_bartel/nbisaria/analysis/%s/%s/%s/%s",
	                     mirna, experiment, library_type, kd_type)
	full_path <- file.path(directory,
	                       "kds_singlebg_multinomial_gradient_pars.txt")

	dG_table <- read.table(full_path, sep=" ", row.names=1, skip=1, header=FALSE)
	colnames(dG_table) <- gsub("^(X)", colnames(dG_table), replace="")
	dG_table
}

Get8merSeedMismatchTypes <- function(mirna) {
	site_8mer <- GetSiteSeq(mirna, "8mer")
	positions_vec <- unlist(strsplit(site_8mer, split=""))
	mm_names_all <- c()
	for (i in 2:7) {
		nuc_i <- positions_vec[9 - i]
		mm_nucs <- setdiff(kDNucs, nuc_i)
		if (nuc_i == "C") {
			mm_names <- c(sprintf("mm%s%s", c("A", "G"), i), sprintf("w%s", i))
		} else if (nuc_i == "A") {
			mm_names <- c(sprintf("mm%s%s", c("C", "T"), i), sprintf("w%s", i))
		} else {
			mm_names <- sprintf("mm%s%s", mm_nucs, i)
		}
		mm_names_all <- c(mm_names_all, mm_names)
	}
	seqs <- sapply(mm_names_all, function(mm_name) {
		GetSiteSeq(mirna, sprintf("8mer-%s", mm_name))
	})
	site_names <- sprintf("8mer-%s",  mm_names_all)
	names(site_names) <- seqs
	site_names
}

GetInternalMismatches <- function(threep_site_name, threep_site_seq, mirna) {
	len_site <- nchar(threep_site_seq)
	positions_vec <- unlist(strsplit(threep_site_seq, split=""))
	mm_names_all <- c()
	mirna_start_pos <- as.integer(unlist(
    strsplit(
			unlist(
				strsplit(
					threep_site_name, split="\\."
				)
			)[1],
			split="-m"
		)
	)[2])

	# Address 3′ nucleotide
	i = 1
	nuc_i <- positions_vec[len_site + 1 - i]
	mm_nucs <- setdiff(kDNucs, nuc_i)
	if (nuc_i == "C") {
		mm_names <- c(sprintf("%smer-m%s.%s-%s", len_site - 1, mirna_start_pos + 1,
		                      mirna_start_pos + len_site - 1, c("A", "G")),
		              sprintf("%sw%s", threep_site_name, i + mirna_start_pos - 1))
	} else if (nuc_i == "A") {
		mm_names <- c(sprintf("%smer-m%s.%s-%s", len_site - 1, mirna_start_pos + 1,
		                      mirna_start_pos + len_site - 1, c("C", "T")),
		              sprintf("%sw%s", threep_site_name, i + mirna_start_pos - 1))
	} else {
		mm_names <- sprintf("%smer-m%s.%s-%s", len_site - 1, mirna_start_pos + 1,
		                      mirna_start_pos + len_site - 1, mm_nucs)
	}
	mm_names_all <- c(mm_names_all, mm_names)

	for (i in 2:(len_site - 1)) {
		nuc_i <- positions_vec[len_site + 1 - i]
		mm_nucs <- setdiff(kDNucs, nuc_i)
		if (nuc_i == "C") {
			mm_names <- c(sprintf("%s-mm%s%s", threep_site_name,
			                      c("A", "G"), i + mirna_start_pos - 1),
			              sprintf("%sw%s", threep_site_name,
			                      i + mirna_start_pos - 1))
		} else if (nuc_i == "A") {
			mm_names <- c(sprintf("%s-mm%s%s", threep_site_name, c("C", "T"),
			                      i + mirna_start_pos - 1),
			              sprintf("%sw%s", threep_site_name,
			                      i + mirna_start_pos - 1))
		} else {
			mm_names <- sprintf("%s-mm%s%s", threep_site_name, mm_nucs,
			                    i + mirna_start_pos - 1)
		}
		mm_names_all <- c(mm_names_all, mm_names)
	}
	# Address 5′ nucleotide
	i = len_site
	nuc_i <- positions_vec[len_site + 1 - i]
	mm_nucs <- setdiff(kDNucs, nuc_i)
	if (nuc_i == "C") {
		mm_names <- c(sprintf("%s-%smer-m%s.%s", c("A", "G"), len_site - 1, mirna_start_pos,
		                      mirna_start_pos + len_site - 2),
		              sprintf("%sw%s", threep_site_name, i + mirna_start_pos - 1))
	} else if (nuc_i == "A") {
		mm_names <- c(sprintf("%s-%smer-m%s.%s", c("C", "T"), len_site - 1, mirna_start_pos,
		                      mirna_start_pos + len_site - 2),
		              sprintf("%sw%s", threep_site_name, i + mirna_start_pos - 1))
	} else {
		mm_names <- sprintf("%s-%smer-m%s.%s", mm_nucs, len_site - 1, mirna_start_pos,
		                      mirna_start_pos + len_site - 2)
	}
	mm_names_all <- c(mm_names_all, mm_names)
	# Get the sequences of these sites using the python function.
	seqs <- sapply(mm_names_all, function(mm_name) {
		GetSiteSeq(mirna, mm_name)
	})
	# Use these sequences as names of the site vector to function as a dictionary.
	names(mm_names_all) <- seqs
	# Return the named vector.
	mm_names_all
}

GetInternalBulges <- function(threep_site_name, threep_site_seq, mirna) {
	len_site <- nchar(threep_site_seq)
	positions_vec <- unlist(strsplit(threep_site_seq, split=""))
	mm_names_all <- c()
	mirna_start_pos <- as.integer(unlist(
    strsplit(
			unlist(
				strsplit(
					threep_site_name, split="\\."
				)
			)[1],
			split="-m"
		)
	)[2])
	bulge_names_all <- c()
	for (i in 1:(len_site)) {
		for (nuc in kDNucs) {
			bulge_name <- sprintf("%sb%s%s", threep_site_name, nuc, i + mirna_start_pos - 1)
			bulge_names_all <- c(bulge_names_all, bulge_name)
		}
	}
	seqs <- sapply(bulge_names_all, function(mm_name) {
		GetSiteSeq(mirna, mm_name)
	})
	dup_bool <- duplicated(seqs)
	dup_seqs <- seqs[dup_bool]
	uniq_seqs <- seqs[!dup_bool]
	bulgepos.names <- uniq_seqs
	pos_list <- vector("list", length(uniq_seqs))
	nuc_list <- vector("list", length(uniq_seqs))
	names(pos_list) <- uniq_seqs
	names(nuc_list) <- uniq_seqs
	for (i in 1:length(names(uniq_seqs))) {
		name <- names(uniq_seqs)[i]
		name_num <- unlist(strsplit(name, split="b"))[2]
		num <- as.integer(substr(name_num, 2, nchar(name_num)))
		nuc <- substr(name_num, 1, 1)
		pos_list[[uniq_seqs[i]]] <- c(num)
		nuc_list[[uniq_seqs[i]]] <- nuc
	}
	for (i in 1:length(dup_seqs)) {
		dup_seq <- dup_seqs[i]
		dup_name <- names(dup_seqs)[i]
		dup_name_num <- unlist(strsplit(dup_name, split="b"))[2]
		dup_num <- as.integer(substr(dup_name_num, 2, nchar(dup_name_num)))
		pos_list[[dup_seq]] <- c(pos_list[[dup_seq]], dup_num)
	}
	names_final <- vector(length=length(uniq_seqs))
	for (i in 1:length(names_final)) {
		nuc <- nuc_list[[uniq_seqs[i]]]
		pos <- pos_list[[uniq_seqs[i]]]
		len_bulge <- length(pos)
		if (len_bulge > 1) {
			range <- c(min(pos), max(pos))
			range_string <- sprintf("(%s.%s)", range[1], range[2])
		} else {
			range_string <- sprintf("%s", pos)
		}
		name <- sprintf("%sb%s%s", threep_site_name, nuc, range_string)
		names_final[i] <- name
	}
	names(uniq_seqs) <- names_final
	inds_keep <- grep(sprintf("%s$", mirna_start_pos),
	                  names(uniq_seqs), perl=TRUE, invert=TRUE)
	seqs_final <- uniq_seqs[inds_keep]

	names_final <- names(seqs_final)
	names(names_final) <- seqs_final
	names_final
}

GetInternalDeletions <- function(threep_site_name, threep_site_seq, mirna) {
	len_site <- nchar(threep_site_seq)
	positions_vec <- unlist(strsplit(threep_site_seq, split=""))
	mm_names_all <- c()
	mirna_start_pos <- as.integer(unlist(
    strsplit(
			unlist(
				strsplit(
					threep_site_name, split="\\."
				)
			)[1],
			split="-m"
		)
	)[2])
	del_names_all <- c()
	for (i in 1:(len_site)) {
		del_name <- sprintf("%sd%s", threep_site_name, i + mirna_start_pos - 1)
		del_names_all <- c(del_names_all, del_name)
	}
	seqs <- sapply(del_names_all, function(mm_name) {
		GetSiteSeq(mirna, mm_name)
	})
	dup_bool <- duplicated(seqs)
	dup_seqs <- seqs[dup_bool]
	uniq_seqs <- seqs[!dup_bool]
	delpos.names <- uniq_seqs
	pos_list <- vector("list", length(uniq_seqs))
	# nuc_list <- vector("list", length(uniq_seqs))
	names(pos_list) <- uniq_seqs
	# names(nuc_list) <- uniq_seqs
	for (i in 1:length(names(uniq_seqs))) {
		name <- names(uniq_seqs)[i]
		num <- as.integer(unlist(strsplit(name, split="d"))[2])
		# num <- as.integer(substr(name_num, 1, nchar(name_num)))
		# nuc <- substr(name_num, 1, 1)
		pos_list[[uniq_seqs[i]]] <- c(num)
		# nuc_list[[uniq_seqs[i]]] <- nuc
	}
	# break
	if (length(dup_seqs) != 0) {
		for (i in 1:length(dup_seqs)) {
			dup_seq <- dup_seqs[i]
			dup_name <- names(dup_seqs)[i]
			dup_num <- as.integer(unlist(strsplit(dup_name, split="d"))[2])
			# dup_num <- as.integer(substr(dup_name_num, 2, nchar(dup_name_num)))
			pos_list[[dup_seq]] <- c(pos_list[[dup_seq]], dup_num)
		}		
	}
	names_final <- vector(length=length(uniq_seqs))
	for (i in 1:length(names_final)) {
		# nuc <- nuc_list[[uniq_seqs[i]]]
		pos <- pos_list[[uniq_seqs[i]]]
		len_del <- length(pos)
		if (len_del > 1) {
			range <- c(min(pos), max(pos))
			range_string <- sprintf("(%s.%s)", range[1], range[2])
		} else {
			range_string <- sprintf("%s", pos)
		}
		name <- sprintf("%sd%s", threep_site_name, range_string)
		names_final[i] <- name
	}
	names(uniq_seqs) <- names_final
	inds_keep <- grep(sprintf("%s$", mirna_start_pos),
	                  names(uniq_seqs), perl=TRUE, invert=TRUE)
	seqs_final <- uniq_seqs[inds_keep]

	names_final <- names(seqs_final)
	names(names_final) <- seqs_final
	names_final
}

Get3PrimeSequences <- function(mirna, site_len, bulge=FALSE) {
	mirna_len <- nchar(kMirnaSeqs[mirna])
	mirna_starts <- seq(9, (mirna_len - site_len + 1))
	site_names <- sprintf("%smer-m%s.%s", site_len, mirna_starts,
	                      mirna_starts + site_len - 1)
	site_seqs <- sapply(site_names, GetSiteSeq, mirna=mirna)
	site_names <- names(site_seqs)
	names(site_names) <- site_seqs
	for (ind in seq(length(site_names))) {
		site_mm <- GetInternalMismatches(site_names[ind], names(site_names)[ind],
		                                 mirna)
		site_names <- c(site_names, site_mm)
		if (bulge) {
			site_bulge <- GetInternalBulges(site_names[ind], names(site_names)[ind],
			                                mirna)
			site_names <- c(site_names, site_bulge)
			if (site_len > 5) {
				site_del <- GetInternalDeletions(site_names[ind], names(site_names)[ind],
				                                mirna)				
				site_names <- c(site_names, site_del)
			}
		}
	}
	site_names
}

GetKdProperties <- function(kd_list, mirna, bulge=FALSE) {
	# This function takes the rownames of a kd_list generated by Namita
	# and parses the sequences such that it assigns what the 3' and seed regions
	# are in terms of standard site nomenclature.
	kd_names <- rownames(kd_list)

	# First change all of the "8mer-1.8"-style names to "8mer"
	# NOTE, instead of dashes this dictionary uses underscores so that the dashes
	# can be used to separate the columns further down in this function.
	seed_dict <- kSiteNameConversion[["let-7a"]]
	seed_inds <- which(kd_names %in% names(seed_dict))
	kd_names[seed_inds] <- sprintf("NA-NA-%s", seed_dict[kd_names[seed_inds]])

	# Get the matrix of three_p_site, loop_length, and seed_sequence.
	sequence_matrix <- matrix(unlist(sapply(kd_names, function(name) {
		strsplit(name, split="-")
	})), ncol=3, byrow=TRUE)

	# Convert the values in the third column (the seed mismatch type) to their
	# appropriate name using the Get8merSeedMismatchTypes, that returns a named
	# list where the names are the sequences, and the entries are the site types.
	seedmm_dict <- Get8merSeedMismatchTypes(mirna)
	seedmm_inds <- which(sequence_matrix[, 3] %in% names(seedmm_dict))
	sequence_matrix[seedmm_inds, 3] <- seedmm_dict[sequence_matrix[seedmm_inds, 3]]

	# Get the list of three-prime strings to swap the sequences in column 1.
	# iteratively move through lists of length 4 to 9
	for (site_length in seq(5, 9)) {
		threep_dict <- Get3PrimeSequences(mirna, site_length, bulge=bulge)
		threep_inds <- which(sequence_matrix[, 1] %in% names(threep_dict))
		sequence_matrix[threep_inds, 1] <- threep_dict[sequence_matrix[threep_inds, 1]]		
		check_3p <- sequence_matrix[, 1]
		print(length(grep("mer", check_3p, invert=TRUE)))	
	}
	return(sequence_matrix)
}


