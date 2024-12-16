source("general/general.R")
# source("general/ThrPFunctions_temp.R")
# source("general/ThrPFunctions.R")
options(width = 160)

graphics.off()

library(rjson)

# kathy r squareds:
r_squareds_path <- "/lab/bartel4_ata/kathyl/RNA_Seq/outputs/resubmission/r2s/biochem_r2s.txt"


GetKathyMirnaName <- function(mirna) {
  alt_names <- c(`lsy-6`="lsy6", `let-7a`="let7")
  if (mirna %in% names(alt_names)) return(alt_names[[mirna]])
  else return(gsub("^(.*)R\\-(.*)$", replacement="\\1r\\2", mirna, perl=TRUE))
}

GetSeanMirnaName <- function(mirnas) {
  sapply(mirnas, function(mirna) {
    alt_names <- c(`lsy6`="lsy-6", `let7`="let-7a")
    if (mirna %in% names(alt_names)) return(alt_names[[mirna]])
    else return(gsub("^(.*)r(.*)$", replacement="\\1R-\\2", mirna, perl=TRUE))
  })
}

ProcessKathyPredictionsFile <- function(model, cell_line, mirna_set, mir_columns=TRUE, norm=FALSE,
                                        data=FALSE, rsquared=FALSE) {
  base_path <- "/lab/solexa_bartel/klin/miRNA_models_data_old/model_outputs/biochem/"
  if (cell_line == "HeLa") {
    if (mirna_set == "five") {
      if (model == "biochem") {
        final_path <- "pred_dfs/resids_orf_utr3_background.txt"
      } else if (model == "biochemplus") {
        final_path <- "pred_dfs/resids_all_feats.txt" 
      }
      col_use <- 4
    } else if (mirna_set == "sixteen") {
      final_path_1 <- "hela_transfection/feat1_allgenes_lr003_nodropout_batch50"
      final_path_2 <- "_rebalancekds2k_noaugment_repweight095_mask_w3_netpred/"
      final_path_3 <- sprintf("train_all_with_passenger_%s_preds.txt", model)
      final_path <- paste0(final_path_1, final_path_2, final_path_3)
      col_use <- 3
    }
    path_file <- paste0(base_path, final_path)
  }
  print(path_file)
  # Get the starting file.
  df <- read.table(path_file, sep="\t", stringsAsFactors=FALSE)
  print(head(df))
  if (rsquared) {
    # print(head(df))
    df_global <<- df
    return(cor(as.numeric(df[-1, 6]), as.numeric(df[-1, 7]))^2)
  }
  if (mir_columns) {
    if (norm) col_use <- col_use + 2
    if (data) col_use <- col_use + 1
    # Identify the genes and mirnas.
    genes <- unique(df[2:nrow(df), 1])
    mirnas <- unique(df[2:nrow(df), 2])
    # Pre-allocate the matrix with the gene and column names.
    df_out <- matrix(NA, nrow=length(genes),
                          ncol=length(mirnas), dimnames=list(genes, mirnas))
    # Iterate over the mirnas to make the columns.
    for (mirna in mirnas) {
      inds_mir <- which(df[, 2] == mirna)
      df_out[df[inds_mir, 1], mirna] <- as.numeric(df[inds_mir, col_use])
    }
    # Rename the columns.
    colnames(df_out) <- sapply(colnames(df_out), GetSeanMirnaName)
    df <- df_out
  }
  return(df)
}

ProcessSeanPredictionsFile <- function(model, cell_line, mirna_set, kds) {
  base_path <- "/lab/solexa_bartel/mcgeary/transfections"
  if (new) new_str <- "_new"
  else     new_str <- ""
  path_file <- sprintf(
    "%s/%s/count_tables/pred_%s_%s_%s%s.txt", base_path, cell_line, model,
    mirna_set, kds, new_str 
  )
  df_out <- data.matrix(read.table(path_file, sep="\t", header=TRUE))
  colnames(df_out) <- gsub("\\.", replacement="-", colnames(df_out))
  return(df_out)
}


# ProcessSeanPredictionsFileNew <- function(
#   model, cell_line, mirna_set, kds, tpm="sm", passenger=FALSE, new=FALSE,
#   get_tpm=FALSE
# ) {
#   base_path <- "/lab/solexa_bartel/mcgeary/transfections"
#   if (passenger) passenger_str <- "_passenger"
#   else           passenger_str <- ""
#   if (new)  new_str <- "_new"
#   else      new_str <- ""
#   path_file <- sprintf(
#     "%s/%s/model_predictions_tf/%s_%s_%s_%stpm%s%s.txt", base_path, cell_line, model,
#     kds, mirna_set, tpm, passenger_str, new_str
#   )
#   print(path_file)
#   df <- read.table(path_file, sep="\t", stringsAsFactors=FALSE, header=TRUE)
#   print(head(df))
#   # Identify the genes and mirnas.
#   genes <- unique(df$transcript)
#   mirnas <- unique(df$mir)
#   # Pre-allocate the matrix with the gene and column names.
#   df_out <- matrix(NA, nrow=length(genes),
#                         ncol=length(mirnas), dimnames=list(genes, mirnas))
#   # Iterate over the mirnas to make the columns.
#   for (mirna in mirnas) {
#     inds_mir <- which(df$mir == mirna)
#     if (get_tpm) vals <- df$label_normed
#     else         vals <- df$pred_normed
#     df_out[df[inds_mir, 1], mirna] <- as.numeric(vals)[inds_mir]
#   }
#   # Rename the columns.
#   colnames(df_out) <- sapply(colnames(df_out), GetSeanMirnaName)
#   return(df_out)
# }


GetBatchNormTransfectionData <- function(cell_line, mirna_set) {
  base_path <- "/lab/solexa_bartel/mcgeary/transfections"
  remaining_path <- sprintf("%s/count_tables/logtpm_batchnormalized.txt", cell_line)
  full_path <- sprintf("%s/%s", base_path, remaining_path)

  full_df <- data.matrix(read.table(full_path))
  colnames(full_df) <- gsub("\\.", replacement="-", x=colnames(full_df))
  if (mirna_set == "sixteen") {
    mirnas_use <- c("lsy-6", "miR-1", "miR-124", "miR-137",
                    "miR-139", "miR-143", "miR-144", "miR-153",
                    "miR-155", "miR-182", "miR-199a", "miR-204","miR-205", "miR-216b", "miR-223", "miR-7")
  } else if (mirna_set == "six") {
    mirnas_use <- c("let-7a", "lsy-6", "miR-1", "miR-124", "miR-155", "miR-7")
  } else if (mirna_set == "five") {
    mirnas_use <- c("lsy-6", "miR-1", "miR-124", "miR-155", "miR-7")
  }
  out_df <- full_df[, mirnas_use]
  return(out_df)
}

MeanCenter <- function(df) {
  return(df - rowMeans(df))
}

ConvertJsonToArray <- function(json_object) {
  out <- c(json_object$freeAGO, json_object$log_decay, json_object$feature_coefs)
  names(out) <- c(json_object$TRAIN_MIRS, "b", json_object$FEATURE_LIST[-1])

  return(out)
}

GetModelFeatureInputFileOld <- function(
  mirna, cell_line, kds, passenger) {
    # These are the original predictions, generated by Kathy.
  path_1 <- file.path("/lab/solexa_bartel/klin/miRNA_models_data_old",
                         "model_inputs/biochem/")
  if (kds == "measured") {
    path_2 <- ("measured_kds/")
  } else if (kds == "predicted") {
    path_2a <- paste0("predicted_kds/feat1_allgenes_lr003_nodropout_batch50_",
                     "rebalancekds2k_noaugment_repweight095_mask_w3_netpred/")
    path_2b <- tolower(cell_line)
    path_2 <- paste0(path_2a, path_2b)
  }
  path <- paste0(path_1, path_2)
  file <- sprintf("%s.txt", GetKathyMirnaName(mirna))
  full_path <- paste0(path, file)
  print(full_path)
  feature_df <- read.table(full_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  feature_df <- feature_df[feature_df$log_kd <= 0, ]

}

GetModelFeatureInputFile <- function(
  mirna, cell_line, kds, passenger, min_pairing=FALSE, offset_opt=FALSE,
  supp_l=FALSE, supp_r=FALSE, offset_tol=FALSE, w_pairing=FALSE, w_offset=FALSE,
  w_supp=FALSE, two_bmodes=FALSE) {
    # These are generated by me.
  path <- file.path("/lab/solexa_bartel/mcgeary/transfections",
                    sprintf("%s/new_feature_files/%s_kds/", cell_line, kds))

  if (min_pairing) {
    path <- sprintf("%smin_pairing_%s,", path, min_pairing)
  }
  if (offset_opt) {
    path <- sprintf("%soffset_opt_%s,", path, offset_opt)
  }
  if (supp_l) {
    path <- sprintf("%ssupp_l_%s,", path, supp_l)
  }
  if (supp_r) {
    path <- sprintf("%ssupp_r_%s,", path, supp_r)
  }
  if (min_pairing) {
    path <- sprintf("%smin_pairing_%s,", path, min_pairing)
  }
  if (min_pairing) {
    path <- sprintf("%smin_pairing_%s,", path, min_pairing)
  }

  path <- gsub("^(.*),$", replacement="\\1/", path)
  file <- sprintf("%s.txt", GetKathyMirnaName(mirna))
  full_path <- paste0(path, file)
  print(full_path)
  feature_df <- read.table(full_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  feature_df <- feature_df[feature_df$log_kd <= 0, ]
}

GetModelPathString <- function (min_pairing=FALSE, offset_opt=FALSE,
  supp_l=FALSE, supp_r=FALSE, offset_tol=FALSE, w_pairing=FALSE, w_offset=FALSE,
  w_supp=FALSE, two_bmodes=FALSE, oligoG=FALSE, model_coefs=FALSE
) {
    mod_string <- ""
  if (min_pairing) {
    mod_string <- sprintf("%smin_pairing_%s,", mod_string, min_pairing)
  }
  if (offset_opt) {
    mod_string <- sprintf("%soffset_opt_%s,", mod_string, offset_opt)
  }
  if (supp_l) {
    mod_string <- sprintf("%ssupp_l_%s,", mod_string, supp_l)
  }
  if (supp_r) {
    mod_string <- sprintf("%ssupp_r_%s,", mod_string, supp_r)
  }
  if (offset_tol) {
    mod_string <- sprintf("%soffset_tol_%s,", mod_string, offset_tol)
  }
  if (w_pairing) {
    mod_string <- sprintf("%sw_pairing_%s,", mod_string, w_pairing)
  }
  if (w_offset) {
    mod_string <- sprintf("%sw_offset_%s,", mod_string, w_offset)
  }
  if (w_supp) {
    mod_string <- sprintf("%sw_supp_%s,", mod_string, w_supp)
  }
  if (two_bmodes) {
    mod_string <- sprintf("%stwo_bmodes,", mod_string)
  }
  if (oligoG) {
    mod_string <- sprintf("%soligoG,", mod_string)
  }
  if (class(model_coefs) == "character") {
    mod_string <- sprintf("%s%s", mod_string, model_coefs)
  }
  mod_string <- gsub("^(.*),$", replacement="\\1", mod_string)
  return(mod_string)
} 

GetModelFit <- function(
  cell_line, model, kds, mirnas, v_tpms="sm", pars_only=TRUE, passenger=FALSE,
  min_pairing=FALSE, offset_opt=FALSE, supp_l=FALSE, supp_r=FALSE,
  offset_tol=FALSE, w_pairing=FALSE, w_offset=FALSE, w_supp=FALSE,
  two_bmodes=FALSE, model_coefs=FALSE
) {
  path1 <- sprintf("/lab/solexa_bartel/mcgeary/transfections/%s/", cell_line)
  path2 <- sprintf("model_parameters_tf/%s_%s_%s_%stpm", model, kds, mirnas,
                   v_tpms)
  path <- paste0(path1, path2)
  mod_string <- SubfunctionCall(GetModelPathString)
  if (mod_string != "") {
    path <- paste0(path, "_", mod_string)
  }
  path <- gsub("^(.*),$", replacement="\\1/", path)
  full_path <- paste0(path, ".txt")
  json_object <- fromJSON(file=full_path)
  out <- c(json_object$freeAGO, json_object$log_decay, json_object$feature_coefs,
           json_object$final_loss, json_object$r2)
  names(out) <- c(json_object$TRAIN_MIRS, "log_decay",
                  json_object$FEATURE_LIST[-1], "final_loss", "r2")
  if (pars_only) {
    out <- out[1:(length(out) - 2)]
  }
  return(out)
}

GetModelPredictionCor <- function(
  cell_line, model, kds, mirnas, v_tpms="sm", passenger=FALSE, min_pairing=FALSE,
  offset_opt=FALSE, supp_l=FALSE, supp_r=FALSE, offset_tol=FALSE,
  w_pairing=FALSE, w_offset=FALSE, w_supp=FALSE, two_bmodes=FALSE, model_coefs=NULL
) {
  path1 <- sprintf("/lab/solexa_bartel/mcgeary/transfections/%s/", cell_line)
  path2 <- sprintf("model_predictions_tf/%s_%s_%s_%stpm", model, kds, mirnas,
                   v_tpms)
  path <- paste0(path1, path2)
  mod_string <- SubfunctionCall(GetModelPathString)
  if (mod_string != "") {
    print(mod_string)
    path <- paste0(path, "_", mod_string)
  }
  full_path <- paste0(path, ".txt")
  file <- read.table(full_path, sep="\t", header=TRUE)
  # Calculate r2 for the model predictions.
  r2 <- cor(file$pred_normed, file$label_normed)^2
  return(r2)
}

GetModelPredictionData <- function(
  cell_line, model, kds, mirnas, v_tpms="sm", passenger=FALSE, min_pairing=FALSE,
  offset_opt=FALSE, supp_l=FALSE, supp_r=FALSE, offset_tol=FALSE,
  w_pairing=FALSE, w_offset=FALSE, w_supp=FALSE, two_bmodes=FALSE,
  model_coefs=FALSE, norm=FALSE, data=FALSE, mir_columns=FALSE
) {
  path1 <- sprintf("/lab/solexa_bartel/mcgeary/transfections/%s/", cell_line)
  path2 <- sprintf("model_predictions_tf/%s_%s_%s_%stpm", model, kds, mirnas,
                   v_tpms)
  path <- paste0(path1, path2)
  mod_string <- SubfunctionCall(GetModelPathString)
  if (mod_string != "") {
    print(mod_string)
    path <- paste0(path, "_", mod_string)
  }
  full_path <- paste0(path, ".txt")
  df <- read.table(full_path, sep="\t", header=TRUE)
  # Calculate r2 for the model predictions.
  if (mir_columns) {
  # Identify the genes and mirnas.
    genes <- unique(df[2:nrow(df), 1])
    mirnas <- unique(df[2:nrow(df), 2])
    # Pre-allocate the matrix with the gene and column names.
    df_out <- matrix(NA, nrow=length(genes),
                          ncol=length(mirnas), dimnames=list(genes, mirnas))
    # Iterate over the mirnas to make the columns.
    col_use <- 3
    if (norm) col_use <- col_use + 2
    if (data) col_use <- col_use + 1
    for (mirna in mirnas) {
      inds_mir <- which(df[, 2] == mirna)
      df_out[df[inds_mir, 1], mirna] <- as.numeric(df[inds_mir, col_use])
    }
    # Rename the columns.
    colnames(df_out) <- sapply(colnames(df_out), GetSeanMirnaName)
    df <- df_out
    # df <- df[which(rownames(df) != "NM_006772.2"),]
    if (length(mirnas) == 5) {
      df <- df[, c("lsy-6", "miR-1", "miR-124", "miR-155", "miR-7")]
    }
  }
  return(df)
}


GetTargetScanPredictions <- function(
  file_index=1, bounded=TRUE, kathy_generated=FALSE, refit=FALSE
) {
  if (kathy_generated) {
    if (refit)  type_str <- "multisite"
    else        type_str <- "original"
    full_path <- sprintf("/lab/bartel4_ata/kathyl/RNA_Seq/outputs/ts7/no_xval_results/%s_predictions.txt",
                         type_str)
    print(full_path)
  } else {
    full_path <- sprintf("Repression/temp_targetscan_output_%s/predictions.txt",
                         file_index)
  }
  print(full_path)
  df <- read.table(full_path, sep="\t", header=TRUE)
  # Identify the genes and mirnas.
  genes <- unique(df[1:nrow(df), 1])
  mirnas <- unique(df[1:nrow(df), 2])
  # Pre-allocate the matrix with the gene and column names.
  df_out <- matrix(0, nrow=length(genes),
                        ncol=length(mirnas), dimnames=list(genes, mirnas))
  # Iterate over the mirnas to make the columns.
  col_use <- 3
  if (bounded) col_use <- col_use + 1
  for (mirna in mirnas) {
    inds_mir <- which(df[, 2] == mirna)
    df_out[df[inds_mir, 1], mirna] <- as.numeric(df[inds_mir, col_use])
  }
  # Rename the columns.
  colnames(df_out) <- sapply(colnames(df_out), GetSeanMirnaName)
  df <- df_out
  df <- df_out[, colnames(df) != "let-7a"]
  # df <- df[, c("lsy-6", "miR-1", "miR-124", "miR-155", "miR-7")]
  return(df)
}


GetTargetScanPreds <- function(
  cell_line, train_mirnas="sixteen", norm=FALSE, tpms=FALSE,
  rescaled=FALSE, bounded=FALSE, alt=FALSE, min_pairing=FALSE,
  offset_opt=FALSE, supp_l=FALSE, supp_r=FALSE, offset_tol=FALSE,
  w_pairing=FALSE, w_offset=FALSE, w_supp=FALSE, two_bmodes=FALSE, oligoG=FALSE, model_coefs=NULL
) {
  path1 <- sprintf("/lab/solexa_bartel/mcgeary/transfections/%s/", cell_line)
  if (tpms) name <- "tpms"
  else      name <- "predictions"
  path2 <- sprintf("target_scan/%s/%s", name, name)
  if (norm) path2 <- paste0(path2, "_mc")
  path <- paste0(path1, path2)
  mod_string <- SubfunctionCall(GetModelPathString)
  if (mod_string != "") {
    print(mod_string)
    path <- paste0(path, "_", mod_string)
  }
  path <- paste0(path, "_", train_mirnas)
  if (alt) {
    path <- paste0(path, "_alt")    
  }
  if (rescaled) {
    path <- paste0(path, "_rescaled")
  }
  if (bounded) {
    path <- paste0(path, "_bounded")
  }
  full_path <- paste0(path, ".txt")
  print(full_path)
  file <- read.table(full_path, sep="\t", header=TRUE, row.names=1)
  return(file)
}


GetTargetScanFeatureFile <- function(cell_line, min_pairing=FALSE,
  offset_opt=FALSE, supp_l=FALSE, supp_r=FALSE, offset_tol=FALSE,
  w_pairing=FALSE, w_offset=FALSE, w_supp=FALSE, two_bmodes=FALSE, oligoG=FALSE,
  model_coefs=NULL
) {
  path1 <- sprintf("/lab/solexa_bartel/mcgeary/transfections/%s/", cell_line)
  path2 <- "target_scan/new_threep_scores/"
  path <- paste0(path1, path2)
  mod_string <- SubfunctionCall(GetModelPathString)
  if (mod_string != "") {
    print(mod_string)
    path <- paste0(path, mod_string, "/")
  }
  full_path <- paste0(path, "features.txt")
  print(full_path)
  file <- read.table(full_path, sep="\t", header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
  return(file)
}


GetTargetScanR2 <- function(
  cell_line, mirnas, train_mirnas, rescaled=FALSE, bounded=FALSE, alt=FALSE, min_pairing=FALSE,
  offset_opt=FALSE, supp_l=FALSE, supp_r=FALSE, offset_tol=FALSE,
  w_pairing=FALSE, w_offset=FALSE, w_supp=FALSE, two_bmodes=FALSE, model_coefs=NULL
) {
  path1 <- sprintf("/lab/solexa_bartel/mcgeary/transfections/%s/", cell_line)
  path2 <- sprintf("target_scan/r2/r2")
  path <- paste0(path1, path2)
  mod_string <- SubfunctionCall(GetModelPathString)
  if (mod_string != "") {
    print(mod_string)
    path <- paste0(path, "_", mod_string)
  }
  path <- paste0(path, "_", train_mirnas)
  if (alt) {
    path <- paste0(path, "_alt")    
  }
  if (rescaled) {
    path <- paste0(path, "_rescaled")
  }
  if (bounded) {
    path <- paste0(path, "_bounded")
  }
  full_path <- paste0(path, ".txt")
  print(full_path)
  file <- read.table(full_path, sep="\t", header=TRUE, row.names=1)[[mirnas]]
  print(file)
  return(file)
}
# pars_orig <- read.table("Repression/kathy_scripts/target_scan/multisite_params.txt", sep="\t", header=TRUE, row.names=1)
# pars_new <- read.table("Repression/new_parameters.txt", sep="\t", header=TRUE, row.names=1)

# pars_new_norm <- read.table("Repression/brand_new_pars_formatted.txt", sep="\t", header=TRUE, row.names=1)


# print(pars_orig)
# print(pars_new)

# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(unlist(pars_orig), unlist(pars_new))

# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(unlist(pars_orig), unlist(pars_new_norm))

# break

# cell_line <- "HeLa"

# mirnas <- "five"


# five <- c("lsy-6", "miR-1", "miR-124", "miR-155", "miR-7")
# preds_biochem <- GetModelPredictionData(cell_line, "biochem", "measured", mirnas, norm=FALSE, mir_columns=TRUE)
# preds_biochemplus <- GetModelPredictionData(cell_line, "biochemplus", "measured", mirnas, norm=FALSE, mir_columns=TRUE)
# # preds_new <- GetModelPredictionData(cell_line, "biochemplus", "measured", mirnas, two_bmodes=TRUE, norm=TRUE, mir_columns=TRUE)
# data_norm_alt <- GetModelPredictionData(cell_line, "biochemplus", "measured", "six", norm=TRUE, data=TRUE, mir_columns=TRUE)[, colnames(preds_new)]
# data_norm <- GetModelPredictionData(cell_line, "biochemplus", "measured", mirnas, norm=TRUE, data=TRUE, mir_columns=TRUE)

# data_norm_alt <- GetModelPredictionData(cell_line, "biochemplus", "predicted", "sixteen", norm=TRUE, data=TRUE, mir_columns=TRUE)
# data_norm_alt <- data_norm_alt[, colnames(data_norm)]


# preds_ts7_1 <- GetTargetScanPredictions(file_index="orig")
# preds_ts7_2 <- GetTargetScanPredictions(file_index="rescale")
# preds_ts7_3 <- GetTargetScanPredictions(file_index="bounded")

# preds_ts7_1_alt <- GetTargetScanPredictions(kathy_generated=TRUE)
# preds_ts7_1_alt <- preds_ts7_1_alt[, colnames(preds_ts7_1)]

# print(head(preds_ts7_1_alt))
# print(head(preds_ts7_1))

# inds_missing_ts7_1 <- setdiff(rownames(data_norm), rownames(preds_ts7_1))
# zeros_mat_1 <- matrix(0, nrow=length(inds_missing_ts7_1), ncol=ncol(preds_ts7_1),
#                     dimnames=list(inds_missing_ts7_1, colnames(preds_ts7_1)))
# preds_ts7_1 <- rbind(preds_ts7_1, zeros_mat_1)
# preds_ts7_1 <- preds_ts7_1[rownames(data_norm), ] 

# inds_missing_ts7_2 <- setdiff(rownames(data_norm), rownames(preds_ts7_2))
# zeros_mat_2 <- matrix(0, nrow=length(inds_missing_ts7_2), ncol=ncol(preds_ts7_2),
#                     dimnames=list(inds_missing_ts7_2, colnames(preds_ts7_2)))
# preds_ts7_2 <- rbind(preds_ts7_2, zeros_mat_2)
# preds_ts7_2 <- preds_ts7_2[rownames(data_norm), ] 

# inds_missing_ts7_3 <- setdiff(rownames(data_norm), rownames(preds_ts7_3))
# zeros_mat_3 <- matrix(0, nrow=length(inds_missing_ts7_3), ncol=ncol(preds_ts7_3),
#                     dimnames=list(inds_missing_ts7_3, colnames(preds_ts7_3)))
# preds_ts7_3 <- rbind(preds_ts7_3, zeros_mat_3)
# preds_ts7_3 <- preds_ts7_3[rownames(data_norm), ] 



# preds_ts7_2 <- GetTargetScanPredictions(file_index=2)
# preds_ts7_2 <- GetTargetScanPredictions(kathy_generated=TRUE)

# preds_ts7_3 <- GetTargetScanPredictions(file_index="rescale")

# inds_missing_ts7_3 <- setdiff(rownames(data_norm), rownames(preds_ts7_3))
# zeros_mat <- matrix(0, nrow=length(inds_missing_ts7_3), ncol=ncol(preds_ts7_3),
#                     dimnames=list(inds_missing_ts7_3, colnames(preds_ts7_3)))
# preds_ts7_3 <- rbind(preds_ts7_3, zeros_mat)
# preds_ts7_3 <- preds_ts7_3[rownames(data_norm), ] 



# preds_ts7_4 <- GetTargetScanPredictions(file_index=4)

# inds_missing_ts7_4 <- setdiff(rownames(data_norm), rownames(preds_ts7_4))
# zeros_mat <- matrix(0, nrow=length(inds_missing_ts7_4), ncol=ncol(preds_ts7_4),
#                     dimnames=list(inds_missing_ts7_4, colnames(preds_ts7_4)))
# preds_ts7_4 <- rbind(preds_ts7_4, zeros_mat)
# preds_ts7_4 <- preds_ts7_4[rownames(data_norm), ] 


# preds_ts7_5 <- GetTargetScanPredictions(kathy_generated=TRUE, refit=TRUE)


# preds_ts7_6 <- GetTargetScanPredictions(file_index=5)

# inds_missing_ts7_6 <- setdiff(rownames(data_norm), rownames(preds_ts7_6))
# zeros_mat <- matrix(0, nrow=length(inds_missing_ts7_6), ncol=ncol(preds_ts7_6),
#                     dimnames=list(inds_missing_ts7_6, colnames(preds_ts7_6)))
# preds_ts7_6 <- rbind(preds_ts7_6, zeros_mat)
# preds_ts7_6 <- preds_ts7_6[rownames(data_norm), ] 


# threep_target_orig <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/new_threep_scores/features.txt",
#                                  sep="\t", header=TRUE, row.names=NULL, stringsAsFactors=FALSE)

# threep_target_twobmodes <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/new_threep_scores/two_bmodes/features.txt",
#                                  sep="\t", header=TRUE, row.names=NULL, stringsAsFactors=FALSE)

# threep_target_minpairing3 <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/new_threep_scores/min_pairing_3/features.txt",
#                                  sep="\t", header=TRUE, row.names=NULL, stringsAsFactors=FALSE)



# inds_use <- which(threep_target_orig[, "miRNA.family"] == "mir124")

# print(head(threep_target_orig))
# print(head(threep_target_twobmodes))

# pairings_orig <- threep_target_orig[inds_use, "Threep_pairing"]
# scores_orig <- threep_target_orig[inds_use, "Threep.score"]
# pairings_twobmodes <- threep_target_twobmodes[inds_use, "Threep_pairing"]
# scores_twobmodes <- threep_target_twobmodes[inds_use, "Threep.score"]


# inds_new <- which(pairings_orig != pairings_twobmodes)

# print(cbind(pairings_orig, scores_orig, pairings_twobmodes, scores_twobmodes)[inds_new, ])
# break
# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(threep_target_orig[inds_use, "Threep.score"], threep_target_twobmodes[inds_use, "Threep.score"])

# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(threep_target_orig[inds_use, "Threep.score"], threep_target_minpairing3[inds_use, "Threep.score"])


# break


# pars_orig <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/params/params.txt", sep="\t", header=TRUE, row.names=1)
# pars_rescale <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/params/params_rescaled.txt", sep="\t", header=TRUE, row.names=1)


# preds_orig <- MeanCenter(read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions/predictions.txt", sep="\t", header=TRUE, row.names=1))
# preds_rescale <- MeanCenter(read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions/predictions_rescale.txt", sep="\t", header=TRUE, row.names=1))
# # preds_alt <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions/predictions_alt.txt", sep="\t", header=TRUE, row.names=1)
# # preds_alt_rescale <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions/predictions_alt_rescale.txt", sep="\t", header=TRUE, row.names=1)
# preds_bound <- MeanCenter(read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions/predictions_bounded.txt", sep="\t", header=TRUE, row.names=1))
# # preds_bound_rescale <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions/predictions_bounded_rescale.txt", sep="\t", header=TRUE, row.names=1)

# tpms_orig <- MeanCenter(read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/tpms/tpms.txt", sep="\t", header=TRUE, row.names=1))
# # colnames(preds_orig) <- GetSeanMirnaName(colnames(preds_orig))
# # preds_orig <- preds_orig[, colnames(preds_ts7_1)]
# # print(head(preds_orig))
# # print(head(preds_ts7_1))
# # break
# # print(head(preds_rescale))

# preds_orig_sep <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions_separated/predictions_separated.txt", sep="\t", header=TRUE, row.names=NULL)
# preds_rescale_sep <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions_separated/predictions_separated_rescaled.txt", sep="\t", header=TRUE, row.names=NULL)

# print(head(preds_orig_sep))
# print(head(preds_rescale_sep))


# dev.new(xpos=0, ypos=0, height=3, width=3)
# par(mar=c(2, 2, 2, 1))

# col_sites <- c("purple", "red", "blue", "cyan")

# site_vals <- unique(preds_orig_sep[, "Stype"])
# inds_6mer <- which(preds_orig_sep[, "Stype"] == site_vals[1])
# inds_7merA1 <- which(preds_orig_sep[, "Stype"] == site_vals[2])
# inds_7merm8 <- which(preds_orig_sep[, "Stype"] == site_vals[3])
# inds_8mer <- which(preds_orig_sep[, "Stype"] == site_vals[4])
# cols <- rep("black", nrow(preds_orig_sep))
# cols[inds_6mer] <- ConvertRColortoRGB(kSiteColors["6mer"], alpha=0.5)
# cols[inds_7merA1] <- ConvertRColortoRGB(kSiteColors["7mer-A1"], alpha=0.5)
# cols[inds_7merm8] <- ConvertRColortoRGB(kSiteColors["7mer-m8"], alpha=0.5)
# cols[inds_8mer] <- ConvertRColortoRGB(kSiteColors["8mer"], alpha=0.5)

# x <- rowSums(preds_orig_sep[, 3:ncol(preds_orig_sep)])
# y <- rowSums(preds_rescale_sep[, 3:ncol(preds_rescale_sep)])
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="Summed score")


# dev.new(xpos=300, ypos=0, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "TA"]
# y <- preds_rescale_sep[, "TA"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="Target site abundance")
# dev.new(xpos=600, ypos=0, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "SPS"]
# y <-  preds_rescale_sep[, "SPS"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="Seep pairing stability")

# dev.new(xpos=900, ypos=0, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "SA"]
# y <- preds_rescale_sep[, "SA"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="Site Accessibility")
# dev.new(xpos=1200, ypos=0, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- rowSums(preds_orig_sep[, c("siRNA_1A", "siRNA_1C")])
# y <- rowSums(preds_rescale_sep[, c("siRNA_1A", "siRNA_1C")])
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="miRNA pos. 1")
# dev.new(xpos=1500, ypos=0, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- rowSums(preds_orig_sep[, c("siRNA_8C", "siRNA_8A", "siRNA_8G")])
# y <- rowSums(preds_rescale_sep[, c("siRNA_8C", "siRNA_8A", "siRNA_8G")])
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="miRNA pos. 8")
# dev.new(xpos=0, ypos=300, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "Local_AU_score"]
# y <- preds_rescale_sep[, "Local_AU_score"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="local AU")
# dev.new(xpos=300, ypos=300, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "Threep_score"]
# y <- preds_rescale_sep[, "Threep_score"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="ThreeP score")
# dev.new(xpos=600, ypos=300, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "Min_dist_score"]
# y <- preds_rescale_sep[, "Min_dist_score"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="Min dist. to 3' UTR end")
# dev.new(xpos=900, ypos=300, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "UTR_length_score"]
# y <- preds_rescale_sep[, "UTR_length_score"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="UTR length")
# dev.new(xpos=1200, ypos=300, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "ORF_length"]
# y <- preds_rescale_sep[, "ORF_length"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="ORF length")
# dev.new(xpos=1500, ypos=300, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "ORF_8mers"]
# y <- preds_rescale_sep[, "ORF_8mers"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="ORF 8mer")
# dev.new(xpos=0, ypos=600, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "Off6m_score"]
# y <- preds_rescale_sep[, "Off6m_score"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="Offset 6mer")
# dev.new(xpos=300, ypos=600, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "PCT"]
# y <- preds_rescale_sep[, "PCT"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="PCT")
# dev.new(xpos=600, ypos=600, height=3, width=3)
# par(mar=c(2, 2, 2, 1))
# x <- preds_orig_sep[, "Stype"]
# y <- preds_rescale_sep[, "Stype"]
# lims <- c(min(x, y), max(x, y))
# plot(x, y, xlim=lims, ylim=lims, col=cols)
# segments(x0=lims[1], y0=lims[1], x1=lims[2], y1=lims[2], lty=2)
# segments(x0=lims[1], y0=0, x1=lims[2], col="gray50")
# segments(x0=0, y0=lims[1], y1=lims[2], col="gray50")
# title(main="Stype")

# break
# preds_orig_use <- preds_orig[, c("lsy6", "mir1", "mir124", "mir155", "mir7")]
# preds_rescale_use <- preds_rescale[, c("lsy6", "mir1", "mir124", "mir155", "mir7")]
# # preds_alt_use <- preds_alt[, c("lsy6", "mir1", "mir124", "mir155", "mir7")]
# # preds_alt_rescale_use <- preds_alt_rescale[, c("lsy6", "mir1", "mir124", "mir155", "mir7")]
# preds_bound_use <- preds_bound[, c("lsy6", "mir1", "mir124", "mir155", "mir7")]
# # preds_bound_rescale_use <- preds_bound_rescale[, c("lsy6", "mir1", "mir124", "mir155", "mir7")]

# tpms_use <- tpms_orig[, c("lsy6", "mir1", "mir124", "mir155", "mir7")]

# colnames(preds_orig_use) <- GetSeanMirnaName(colnames(preds_orig_use))
# colnames(preds_rescale_use) <- GetSeanMirnaName(colnames(preds_rescale_use))
# colnames(preds_bound_use) <- GetSeanMirnaName(colnames(preds_bound_use))

# colnames(tpms_use) <- GetSeanMirnaName(colnames(tpms_use))

# preds_orig_use <- preds_orig_use[, five]
# preds_rescale_use <- preds_rescale_use[, five]
# preds_bound_use <- preds_bound_use[, five]

# tpms_use <- tpms_use[, five]


# preds_mc_orig_use_alt <- MeanCenter(preds_ts7_1)
# print(head(preds_mc_orig_use_alt))

# print(cor(unlist(preds_orig_use), unlist(tpms_use))^2)
# print(cor(unlist(preds_rescale_use), unlist(tpms_use))^2)
# # print(cor(unlist(preds_alt_use), unlist(tpms_use))^2)
# # print(cor(unlist(preds_alt_rescale_use), unlist(tpms_use))^2)
# print(cor(unlist(preds_bound_use), unlist(tpms_use))^2)
# # print(cor(unlist(preds_bound_rescale_use), unlist(tpms_use))^2)

# pred_1_norm <- MeanCenter(preds_ts7_1)[, five]
# pred_2_norm <- MeanCenter(preds_ts7_2)[, five]
# pred_3_norm <- MeanCenter(preds_ts7_3)[, five]



# r_squared_1 <- cor(c(pred_1_norm), c(data_norm_alt), use="pairwise.complete.obs")^2
# r_squared_2 <- cor(c(pred_2_norm), c(data_norm_alt), use="pairwise.complete.obs")^2
# r_squared_3 <- cor(c(pred_3_norm), c(data_norm_alt), use="pairwise.complete.obs")^2

# r_squared_1_alt_tpm <- cor(c(pred_1_norm), unlist(tpms_use), use="pairwise.complete.obs")^2
# r_squared_2_alt_tpm <- cor(c(pred_2_norm), unlist(tpms_use), use="pairwise.complete.obs")^2
# r_squared_3_alt_tpm <- cor(c(pred_3_norm), unlist(tpms_use), use="pairwise.complete.obs")^2



# # r_squared_6 <- cor(c(MeanCenter(preds_ts7_6)[, five]), c(data_norm_alt), use="pairwise.complete.obs")^2


# print(r_squared_1)
# print(r_squared_2)
# print(r_squared_3)

# print(r_squared_1_alt_tpm)
# print(r_squared_2_alt_tpm)
# print(r_squared_3_alt_tpm)

# print(cor(unlist(preds_orig_use), c(pred_1_norm))^2)
# print(cor(unlist(preds_orig_use), unlist(preds_rescale_use))^2)

# print(cor(c(data_norm_alt), unlist(tpms_use))^2)
# # print(r_squared_6)
# break

# r_squareds <- c()

# inds_use_1 <- intersect(rownames(preds_ts7_1), rownames(data_norm))
# r_squared_new <- cor(c(MeanCenter(preds_ts7_1[, five])), c(data_norm), use="pairwise.complete.obs")^2
# r_squared_new_alt <- cor(c(MeanCenter(preds_ts7_1)[, five]), c(data_norm_alt), use="pairwise.complete.obs")^2

# print(r_squared_new)
# print(r_squared_new_alt)
# r_squareds <- c(r_squareds, r_squared_new_alt)
# dev.new(xpos=20, ypos=20, height=3.5, width=3.5)
# plot(c(-1*preds_ts7_1[, five]), c(data_norm[, five])/log(2), type="p", xlim=c(0, 2.5), ylim=c(-6, 4))

# r_squared_new <- cor(c(MeanCenter(preds_ts7_2[, five])), c(data_norm), use="pairwise.complete.obs")^2
# r_squared_new_alt <- cor(c(MeanCenter(preds_ts7_2)[, five]), c(data_norm_alt), use="pairwise.complete.obs")^2
# print(r_squared_new)
# print(r_squared_new_alt)
# # r_squareds <- c(r_squareds, r_squared_new)
# dev.new(xpos=20, ypos=320, height=3.5, width=3.5)
# plot(c(-1*preds_ts7_2[, five]), c(data_norm[, five])/log(2), type="p", xlim=c(0, 2.5), ylim=c(-6, 4))

# # inds_use_3 <- intersect(rownames(preds_ts7_3), rownames(data_norm))
# r_squared_new <- cor(c(MeanCenter(preds_ts7_3[, five])), c(data_norm), use="pairwise.complete.obs")^2
# r_squared_new_alt <- cor(c(MeanCenter(preds_ts7_3)[, five]), c(data_norm_alt), use="pairwise.complete.obs")^2
# r_squareds <- c(r_squareds, r_squared_new_alt)
# print(r_squared_new)
# print(r_squared_new_alt)
# dev.new(xpos=320, ypos=20, height=3.5, width=3.5)
# plot(c(-1*preds_ts7_3[, five]), c(data_norm[, five])/log(2), type="p", xlim=c(0, 2.5), ylim=c(-6, 4))

# # inds_use_4 <- intersect(rownames(preds_ts7_4), rownames(data_norm))
# r_squared_new <- cor(c(MeanCenter(preds_ts7_4[, five])), c(data_norm), use="pairwise.complete.obs")^2
# r_squared_new_alt <- cor(c(MeanCenter(preds_ts7_4)[, five]), c(data_norm_alt), use="pairwise.complete.obs")^2
# print(r_squared_new)
# print(r_squared_new_alt)
# dev.new(xpos=320, ypos=320, height=3.5, width=3.5)
# plot(c(-1*preds_ts7_4[, five]), c(data_norm[, five])/log(2), type="p", xlim=c(0, 2.5), ylim=c(-6, 4))


# # inds_use_4 <- intersect(rownames(preds_ts7_4), rownames(data_norm))
# r_squared_new <- cor(c(MeanCenter(preds_biochem)), c(data_norm), use="pairwise.complete.obs")^2
# r_squareds <- c(r_squareds, r_squared_new)
# dev.new(xpos=620, ypos=20, height=3.5, width=3.5)
# plot(c(-1*preds_biochem), c(data_norm)/log(2), type="p", xlim=c(0, 2.5), ylim=c(-6, 4))

# r_squared_new <- cor(c(MeanCenter(preds_biochemplus)), c(data_norm), use="pairwise.complete.obs")^2
# r_squareds <- c(r_squareds, r_squared_new)
# dev.new(xpos=620, ypos=320, height=3.5, width=3.5)
# plot(c(-1*preds_biochemplus), c(data_norm)/log(2), type="p", xlim=c(0, 2.5), ylim=c(-6, 4))

# dev.new(xpos=920, ypos=20, height=3.5, width=3.5)

# xlim <- c(0, 5)
# ylim <- c(0, 0.4)
# par(mar=c(1, 2, 1, 1))
# plot(1, type="n", xlim=xlim, ylim=ylim, axes=FALSE, ann=FALSE)
# axis(side=2, at=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4))
# rect(xleft=0:3, ybottom=0, xright=1:4, ytop=r_squareds, col=c("gray10", "gray30", "goldenrod", "blue"), lwd=0)

# names(r_squareds) <- c("TS7", "TS7_refit", "biochem", "biochemplus")

# # barplot(r_squareds)


# break
# print(diag(cor(preds_orig, data_norm)^2))
# print(diag(cor(preds_new, data_norm)^2))
# print(diag(cor(preds_new, data_norm)^2) - diag(cor(preds_orig, data_norm)^2))

# print(diag(cor(preds_orig, data_norm_alt)^2))
# print(diag(cor(preds_new, data_norm_alt)^2))
# print(diag(cor(preds_new, data_norm_alt)^2) - diag(cor(preds_orig, data_norm_alt)^2))


# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(1:ncol(preds_orig), diag(cor(preds_orig, data_norm)^2), type="o", col="blue")
# points(1:ncol(preds_orig), diag(cor(preds_new, data_norm)^2), type="o", col="red")



# break

# dif_old_new <- (preds_orig[, "lsy-6"] -  preds_new[, "lsy-6"])^2

# rank <- order(dif_old_new, decreasing=TRUE)
# print(head(rank))
# print(preds_orig)
# rank_max <- 1000
# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(preds_orig[rank[1:rank_max], "lsy-6"], preds_new[rank[1:rank_max], "lsy-6"])

# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(preds_orig[rank[1:rank_max], "lsy-6"], data_norm[rank[1:rank_max], "lsy-6"])
# text(x=-0.4, y=0.4, labels=cor(preds_orig[rank[1:rank_max], "lsy-6"], data_norm[rank[1:rank_max], "lsy-6"])^2)

# dev.new(xpos=1020, ypos=20, height=5, width=5)
# plot(preds_new[rank[1:rank_max], "lsy-6"], data_norm[rank[1:rank_max], "lsy-6"])
# text(x=-0.4, y=0.4, labels=cor(preds_new[rank[1:rank_max], "lsy-6"], data_norm[rank[1:rank_max], "lsy-6"])^2)


# temp_five_m_biochem <- GetModelFit(cell_line, "biochem", "measured", "five")
# temp_five_p_biochem <- GetModelFit(cell_line, "biochem", "predicted", "five")

# temp_five_m_biochem_data <- GetModelPredictionData(cell_line, "biochem", "measured", "five")
# temp_five_p_biochem_data <- GetModelPredictionData(cell_line, "biochem", "predicted", "five")

# temp_five_m_biochemplus_data <- GetModelPredictionData(cell_line, "biochemplus", "measured", "five")
# temp_five_p_biochemplus_data <- GetModelPredictionData(cell_line, "biochemplus", "predicted", "five")


# dev.new(xpos=20, ypos=20, height=5, width=5)
# # plot(temp_five_m_biochem_data$pred_normed, temp_five_m_biochem_data$label_normed)
# # text(-0.8, 1.2, cor(temp_five_m_biochem_data$pred_normed, temp_five_m_biochem_data$label_normed)^2, adj=c(0, 0.5))
# inds <- which(temp_five_m_biochem_data$mir == "lsy6")
# plot(temp_five_m_biochem_data$pred_normed[inds], temp_five_p_biochem_data$pred_normed[inds], xlim=c(-1, 1), ylim=c(-1, 1))


# dev.new(xpos=520, ypos=20, height=5, width=5)
# # plot(temp_five_p_biochem_data$pred_normed, temp_five_p_biochem_data$label_normed)
# # text(-0.8, 1.2, cor(temp_five_p_biochem_data$pred_normed, temp_five_p_biochem_data$label_normed)^2, adj=c(0, 0.5))
# inds <- which(temp_five_m_biochem_data$mir == "mir1")
# plot(temp_five_m_biochem_data$pred_normed[inds], temp_five_p_biochem_data$pred_normed[inds], xlim=c(-1, 1), ylim=c(-1, 1))


# dev.new(xpos=1020, ypos=20, height=5, width=5)
# # plot(temp_five_m_biochemplus_data$pred_normed, temp_five_m_biochemplus_data$label_normed)
# # text(-0.8, 1.2, cor(temp_five_m_biochemplus_data$pred_normed, temp_five_m_biochemplus_data$label_normed)^2, adj=c(0, 0.5))
# inds <- which(temp_five_m_biochem_data$mir == "mir124")
# plot(temp_five_m_biochem_data$pred_normed[inds], temp_five_p_biochem_data$pred_normed[inds], xlim=c(-1, 1), ylim=c(-1, 1))



# dev.new(xpos=20, ypos=520, height=5, width=5)
# # plot(temp_five_p_biochemplus_data$pred_normed, temp_five_p_biochemplus_data$label_normed)
# # text(-0.8, 1.2, cor(temp_five_p_biochemplus_data$pred_normed, temp_five_p_biochemplus_data$label_normed)^2, adj=c(0, 0.5))
# inds <- which(temp_five_m_biochem_data$mir == "mir155")
# plot(temp_five_m_biochem_data$pred_normed[inds], temp_five_p_biochem_data$pred_normed[inds], xlim=c(-1, 1), ylim=c(-1, 1))



# dev.new(xpos=520, ypos=520, height=5, width=5)
# inds <- which(temp_five_m_biochem_data$mir == "mir7")
# plot(temp_five_m_biochem_data$pred_normed[inds], temp_five_p_biochem_data$pred_normed[inds], xlim=c(-1, 1), ylim=c(-1, 1))

# dev.new(xpos=1020, ypos=520, height=5, width=5)
# inds <- which(temp_five_m_biochem_data$mir == "mir1")

# plot(temp_five_m_biochem_data$pred_normed[inds], temp_five_p_biochem_data$pred_normed[inds], xlim=c(-1, 1), ylim=c(-1, 1))


# break

# temp_five_m <- GetModelFit(cell_line, "biochemplus", "measured", "five")
# temp_five_p <- GetModelFit(cell_line, "biochemplus", "predicted", "five")
# temp_six_m <- GetModelFit(cell_line, "biochemplus", "measured", "six")
# temp_six_p <- GetModelFit(cell_line, "biochemplus", "predicted", "six")
# temp_sixteen_p <- GetModelFit(cell_line, "biochemplus", "predicted", "sixteen")
# temp_sixteen_p <- GetModelFit(cell_line, "biochemplus", "predicted", "seventeen")


# break

# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, two_bmodes=TRUE))

# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, min_pairing=1))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, min_pairing=3))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, min_pairing=4))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, min_pairing=5))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, min_pairing=6))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, min_pairing=7))

# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_opt=1))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_opt=2))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_opt=3))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_opt=4))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_opt=5))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_opt=6))

# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=9))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=10))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=11))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=12))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_r=17))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_r=18))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_r=19))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_r=20))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_r=21))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_r=22))

# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=11, supp_r=14))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=12, supp_r=15))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=14, supp_r=17))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=15, supp_r=18))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=16, supp_r=19))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, supp_l=17, supp_r=20))

# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_tol=0))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_tol=1))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_tol=2))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_tol=3))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_tol=4))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_tol=5))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, offset_tol=6))

# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=0.1))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=0.2))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=0.3))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=0.4))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=0.5))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=0.6))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=0.7))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=0.8))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=0.9))
# # print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_pairing=1))

# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=0.1))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=0.2))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=0.3))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=0.4))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=0.5))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=0.6))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=0.7))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=0.8))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=0.9))
# # print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_offset=1))

# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=0.1))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=0.2))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=0.3))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=0.4))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=0.5))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=0.6))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=0.7))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=0.8))
# print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=0.9))
# # print(GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas, passenger=passenger, w_supp=1))




# break


# features_old <- GetModelFeatureInputFileOld("let-7a", "HeLa", "measured")
# features_new <- GetModelFeatureInputFile("let-7a", "HeLa", "measured")
# features_diff <- GetModelFeatureInputFile("let-7a", "HeLa", "measured", min_pairing=3)



# print(head(features_old))
# print(head(features_new))
# print(head(features_diff))
# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(features_old$Threep, features_new$Threep)

# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(features_old$Threep, features_diff$Threep)

# break

# pred_old <- ProcessSeanPredictionsFileNew("biochemplus", "HeLa", "five", "measured")

# pred_new <- ProcessSeanPredictionsFileNew("biochemplus", "HeLa", "five", "measured", new=TRUE)

# tpms_old <- ProcessSeanPredictionsFileNew("biochemplus", "HeLa", "five", "measured", get_tpm=TRUE, new=TRUE)
# tpms_new <- ProcessSeanPredictionsFileNew("biochemplus", "HeLa", "five", "measured", get_tpm=TRUE)

# print(head(pred_old))
# print(head(pred_new))
# print(head(tpms_old))
# print(head(tpms_new))

# for (i in seq(5)) {
#   cor_old <- cor(pred_old[, i], tpms_old[, i])^2
#   cor_new <- cor(pred_new[, i], tpms_new[, i])^2
#   mirna <- colnames(pred_old)[i]
#   print(sprintf("%s: %s %s", mirna, cor_old, cor_new))
# }

# print(cor(c(pred_old), c(tpms_old))^2)
# print(cor(c(pred_new), c(tpms_new))^2)

# break


# kl_bc <- ProcessKathyPredictionsFile("biochem", "HeLa", "five")

# kl_bcp <- ProcessKathyPredictionsFile("biochemplus", "HeLa", "five")

# sm_bc <- ProcessSeanPredictionsFile("biochem", "HeLa", "five", "measured")
# sm_bc_new <- ProcessSeanPredictionsFile("biochem", "HeLa", "five", "measured", new=TRUE)

# sm_bcp_new <- ProcessSeanPredictionsFile("biochemplus", "HeLa", "five", "measured", new=TRUE)

# sm_trans_tpm <- GetBatchNormTransfectionData("HeLa", "five")


# kl_bcp_pred <- ProcessKathyPredictionsFile("biochemplus", "HeLa", "sixteen")
# sm_bcp_pred_new <- ProcessSeanPredictionsFile("biochemplus", "HeLa", "sixteen", "predicted", new=TRUE)


# sm_trans_sixteen_tpm <- GetBatchNormTransfectionData("HeLa", "sixteen")

# sm_bc_six_new <- ProcessSeanPredictionsFile("biochem", "HeLa", "six", "measured", new=TRUE)

# sm_trans_six_tpm <- GetBatchNormTransfectionData("HeLa", "six")


# # One gene is different between them: 


# # print(head(sm_trans_tpm))

# print(head(kl_bc))
# print(head(sm_bc))
# print(head(sm_bc_new))


# print(head(kl_bcp))
# print(head(sm_bcp_new))


# # dev.new(xpos=20, ypos=20, height=5, width=5)
# # plot(c(kl_bc), c(sm_bc_new[rownames(kl_bc), ]))

# # dev.new(xpos=520, ypos=20, height=5, width=5)
# # plot(c(kl_bcp), c(sm_bcp_new[rownames(kl_bc), ]))

# # dev.new(xpos=1020, ypos=20, height=5, width=5)
# # plot(c(kl_bcp_pred), c(sm_bcp_pred_new[rownames(kl_bcp_pred), ]))





# dev.new(xpos=20, ypos=520, height=5, width=5)
# rows_use <- intersect(rownames(sm_bc_new), rownames(sm_trans_tpm))
# plot(c(MeanCenter(sm_trans_tpm[rows_use, ])), c(MeanCenter(sm_bc_new[rows_use, ])))
# print(cor(c(MeanCenter(sm_trans_tpm[rows_use, ])), c(MeanCenter(sm_bc_new[rows_use, ])))^2)

# dev.new(xpos=520, ypos=520, height=5, width=5)
# rows_use <- intersect(rownames(sm_bc_six_new), rownames(sm_trans_six_tpm))
# plot(c(MeanCenter(sm_trans_six_tpm[rows_use, ])), c(MeanCenter(sm_bc_six_new[rows_use, ])))
# print(cor(c(MeanCenter(sm_trans_six_tpm[rows_use, ])), c(MeanCenter(sm_bc_six_new[rows_use, ])))^2)

# # dev.new(xpos=520, ypos=520, height=5, width=5)
# # rows_use <- intersect(rownames(sm_bcp_new), rownames(sm_trans_tpm))
# # plot(c(MeanCenter(sm_trans_tpm[rows_use, ])), c(MeanCenter(sm_bcp_new[rows_use, ])))
# # print(cor(c(MeanCenter(sm_trans_tpm[rows_use, ])), c(MeanCenter(sm_bcp_new[rows_use, ])))^2)

# # dev.new(xpos=1020, ypos=520, height=5, width=5)
# # rows_use <- intersect(rownames(sm_bcp_new), rownames(sm_trans_sixteen_tpm))
# # plot(c(MeanCenter(sm_trans_sixteen_tpm[rows_use, ])), c(MeanCenter(sm_bcp_pred_new[rows_use, ])))
# # print(cor(c(MeanCenter(sm_trans_sixteen_tpm[rows_use, ])), c(MeanCenter(sm_bcp_pred_new[rows_use, ])))^2)




# break
# batch_norm_kathy_HeLa <- read.table("/lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/processed/log_tpm_normed.txt")

# batch_norm_sean_HeLa <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/logtpm_batchnormalized.txt")


# batch_norm_sean_HeLa_use <- batch_norm_sean_HeLa[, c("lsy.6", "miR.1", "miR.124", "miR.155", "miR.7")]

# # pred_sm_miR1 <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/pred_miR-1.txt")
# # pred_sm_miR155 <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/pred_miR-155.txt")
# # pred_sm_miR124 <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/pred_miR-124.txt")
# # pred_sm_lsy6 <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/pred_lsy-6.txt")
# # pred_sm_miR7 <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/pred_miR-7.txt")

# # pred_sm_all <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/pred_predicted.txt", sep="\t", header=TRUE)
                          


# kd_bp <- ProcessKathyPredictionsFile
# # pred_sm_miR1 <- read.table("Repression/output_temp/mir1_preds.txt", sep="\t", header=TRUE)

# # rownames(pred_sm_miR1) <- pred_sm_miR1[, 1]
# # pred_sm_miR1 <- pred_sm_miR1[, -1]


# # Print the top of both dataframes in order to figure out how to format them to
# # be similar.

# rownames(pred_sm_all) <- pred_sm_all[, 1]

# pred_sm_all <- pred_sm_all[, 2:6]

# print(head(pred_sm_all))
# print(head(pred_kl_all))

# # Get the unique gene names and the unique miRNA names from kl dataframe; these
# # will be the row and column names, respectively.
# print(kl_matrix_out[1:5, 1:5])

# pred_sm_all <- pred_sm_all[rownames(kl_matrix_out), ]

# print(pred_sm_all[1:5, 1:5])


# print(head(batch_norm_sean_HeLa_use))


# batch_norm_sean_mean_centered <- batch_norm_sean_HeLa_use - rowMeans(batch_norm_sean_HeLa_use)

# batch_norm_sean_mean_centered <- batch_norm_sean_mean_centered[rownames(pred_sm_all), ]

# pred_sm_all_mean_centered <- pred_sm_all - rowMeans(pred_sm_all)


# print(length(unlist(pred_sm_all_mean_centered)))
# print(length(unlist(batch_norm_sean_mean_centered)))


# x_use <- unlist(pred_sm_all_mean_centered)
# y_use <- unlist(batch_norm_sean_mean_centered)

# print(cor(x_use, y_use)^2)




# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(unlist(pred_sm_all), c(kl_matrix_out))

# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(x_use, y_use)

# # plot(kl_matrix_out[, "mir1"], pred_sm_miR1[rownames(pred_sm_all), 3])


# break
# preds_kl_all_use <- pred_kl_all

# break


# # batch_norm_kathy_HEK <- read.table("/lab/solexa_bartel/klin/miRNA_models_data/transfections/hek293ft/processed/log_tpm_normed.txt")

# # batch_norm_sean_HEK <- read.table("/lab/solexa_bartel/mcgeary/transfections/HEK293FT/count_tables/logtpm_batchnormalized.txt")






# print(dim(batch_norm_kathy))
# print(dim(batch_norm_sean))

# print(head(batch_norm_kathy))
# print(head(batch_norm_sean))

# batch_norm_sean_HeLa <- batch_norm_sean_HeLa[, 1:ncol(batch_norm_kathy_HeLa)]

# batch_norm_sean_HEK <- batch_norm_sean_HEK[rownames(batch_norm_kathy_HEK), 1:ncol(batch_norm_kathy_HEK)]

# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(unlist(batch_norm_sean_HeLa), unlist(batch_norm_kathy_HeLa))

# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(unlist(batch_norm_sean_HEK), unlist(batch_norm_kathy_HEK))


# break








# table_sean <- read.table("/lab/solexa_bartel/mcgeary/transfections/HEK293FT/count_tables/raw.txt", stringsAsFactors=FALSE)
# table_kathy <- read.table("/lab/solexa_bartel/klin/miRNA_models_data/transfections/hek293ft/processed/log_tpm.txt")
# colnames(table_kathy) <- gsub("mir", replacement="miR-", colnames(table_kathy))
# colnames(table_kathy) <- gsub("let7", replacement="let-7a", colnames(table_kathy))
# colnames(table_kathy) <- gsub("lsy6", replacement="lsy-6", colnames(table_kathy))

# colnames(table_sean) <- paste0(table_sean[1, ], "_rep", table_sean[2, ])
# rows_use <- intersect(rownames(table_kathy), rownames(table_sean))


# table_sean <- table_sean[rows_use, colnames(table_kathy)]
# table_kathy <- table_kathy[rows_use, ]

# print(table_sean[1:5, 1:5])
# print(table_kathy[1:5, 1:5])

# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(log(c(as.numeric(unlist(table_sean)))), unlist(table_kathy))
# # print(table_sean[1:5, 1:5])



# table_sean <- read.table("/lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/raw.txt", stringsAsFactors=FALSE)
# table_kathy <- read.table("/lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/processed/log_tpm.txt")
# colnames(table_kathy) <- gsub("mir", replacement="miR-", colnames(table_kathy))
# colnames(table_kathy) <- gsub("let7", replacement="let-7a", colnames(table_kathy))
# colnames(table_kathy) <- gsub("lsy6", replacement="lsy-6", colnames(table_kathy))

# colnames(table_sean) <- paste0(table_sean[1, ], "_rep", table_sean[2, ])
# rows_use <- intersect(rownames(table_kathy), rownames(table_sean))


# table_sean <- table_sean[rows_use, colnames(table_kathy)]
# table_kathy <- table_kathy[rows_use, ]

# print(table_sean[1:5, 1:5])
# print(table_kathy[1:5, 1:5])

# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(log(c(as.numeric(unlist(table_sean)))), unlist(table_kathy))



