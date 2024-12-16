bsub python Repression/train.py \
--tpm_file /lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/processed/log_tpm_normed.txt \
--feature_file /lab/solexa_bartel/klin/miRNA_models_data_old/model_inputs/biochem/predicted_kds/feat1_allgenes_lr003_nodropout_batch50_rebalancekds2k_noaugment_repweight095_mask_w3_netpred/hela/MIR.txt \
--mirseqs Repression/mirnas_sixteen.txt \
--kd_cutoff 0 \
--mode all \
--init_bound \
--passenger \
--extra_feats logSA_diff,Threep_canon,PCT \
--outfile temp_outfile_train_sixteen.txt \
--outparams temp_outfile_params_sixteen.txt \
# --tpm_file /lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/logtpm_batchnormalized.txt \
# --feature_file /lab/solexa_bartel/klin/miRNA_models_data_old/model_inputs/biochem/predicted_kds/MIR.txt \
