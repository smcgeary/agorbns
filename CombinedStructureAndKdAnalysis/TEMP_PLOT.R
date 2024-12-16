plot(check/sum(check)*sum(input_temp)/input_temp,check_model/sum(check_model)*sum(input_temp)/input_temp, col = sapply(flanks_temp,GetColorFunction), pch = 19, xlim = c(0,2), ylim = c(0, 2),ann = FALSE, axes = FALSE)

axis(1, at = seq(0, 2, length = 6), lwd = 2, pos = 0)
axis(2, at = seq(0, 2, length = 6), lwd = 2, pos = 0)
title(xlab = "Model enrichment")
title(ylab = "Data enrichment")
dev.new()

plot(check/sum(check)*sum(input_temp)/input_temp,check_model/sum(check_model)*sum(input_temp)/input_temp, col = sapply(flanks_temp,GetColorFunction), log = 'xy', pch = 19, xlim = c(0.2,2), ylim = c(0.2, 2),ann = FALSE, axes = FALSE)

axis(1, at = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2), lwd = 2, pos = 0.2)
axis(2, at = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2), lwd = 2, pos = 0.2)
title(xlab = "Model enrichment")
title(ylab = "Data enrichment")
