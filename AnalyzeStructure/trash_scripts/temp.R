		mean_flanks_predict <- sapply(unique(flanks_I), function(flank) {
		mean(means_I[temp_inds][which(flanks_I[temp_inds]==flank)])/length(temp_inds)
		})