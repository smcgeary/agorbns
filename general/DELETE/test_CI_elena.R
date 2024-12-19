rm()

library(numDeriv)
# source("general/general.R")
# source("general/GenericFigures.R")
graphics.off()
dev.new(xpos=520, ypos=20, height=10, width=10)

## FUNCTIONS

ExponentialFunction <- function(vec_p, vec_t) {
	vec_p[1]*exp(-vec_p[2]*vec_t)
}

# Made it such that ss=TRUE, gives the sum of squared reseiduals;
# ss=FALSE gives the (unsquared) residual vector.
CostFunction <- function(vec_p, vec_t, vec_data, ss=TRUE) {
	vec_y <- ExponentialFunction(vec_p, vec_t)
	if (ss) sum((vec_y - vec_data)^2)
	else vec_y - vec_data
}

# Variance estimation from online resource I found:
EstVar <- function(vec_p, vec_t, vec_data) {
	# OVERALL FORMULA:
	# sum(residuals) / (n - p) * ( A^t %*% A )^(-1)
	# -----a constant---------   ---p x p matrix---
	# Get the matrix A (partial derivatives of each value n by parameter p)
	A <- jacobian(ExponentialFunction, vec_p, vec_t=vec_t)
	# Need to make sigma_2 * ( A^t %*% A )^(-1)
	AtA_inv <- solve(t(A)%*%A)
	# Need n and p for denominator:
	n <- length(vec_t)
	p <- length(vec_p)
	sum_res <- CostFunction(vec_p, vec_t, vec_data)
	sigma_2 <- sum_res/(n - p)
	sigma_2*AtA_inv
}


# Setup underlying model:
vec_t_model <- seq(0, 10, length.out=100)
vec_p <- c(1.5, 0.5)
vec_m <- ExponentialFunction(vec_p, vec_t_model)


# Simulate data:
vec_t <- c(0, 1, 2, 3, 5, 10)
vec_data <- ExponentialFunction(vec_p, vec_t) + rnorm(length(vec_t), 0, 0.07)

# Get the parameter values back by manual non-linear
# least squares regression (optim).
opt <- optim(c(0, 0), CostFunction, vec_t=vec_t, vec_data=vec_data)
vec_p_optim <- opt$par
vec_m_opt <- ExponentialFunction(vec_p_optim, vec_t_model)


model_nls <- nls(vec_data~a*exp(-b*vec_t), start=c(a=1, b=1))
vec_se_optim <- sqrt(diag(EstVar(vec_p_optim, vec_t, vec_data)))

# mat_Var_optim <- EstVar(vec_p_optim, vec_t, vec_data)
# mat_Var_nls <- vcov(model_nls)



print("from optim")
print(vec_p_optim)
print("standard error")
print(vec_se_optim)
print("summary of model")
print(summary(model_nls))

coefs <- coef(summary(model_nls))



vec_se_nls <- coefs[, 2]
vec_p_nls <- coefs[, 1]

CI <- vec_se_nls
par(mfrow=c(2, 2))

plot(vec_t_model, vec_m, type="l")
points(vec_t, vec_data)
lines(vec_t_model, vec_m_opt, type="l", lty=2)

tdist <- qt(0.975, length(vec_t) - length(vec_p_nls))

CIS_manual <- confint(model_nls)

plot(vec_p, vec_p_nls, xlim=c(0, 3), ylim=c(0, 3))
points(vec_p, vec_p_optim, xlim=c(0, 3), ylim=c(0, 3))

arrows(vec_p, vec_p_optim + tdist*vec_se_optim,
       vec_p, vec_p_optim - tdist*vec_se_optim,
       length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3)
arrows(vec_p, vec_p_nls + tdist*vec_se_nls,
       vec_p, vec_p_nls - tdist*vec_se_nls,
       length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3)

arrows(vec_p, CIS_manual[, 1],
       vec_p, CIS_manual[, 2],
       length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3,
       col="red", lty=2)




abline(0, 1, lty=2)
plot(vec_se_optim, vec_se_nls, xlim=c(0, 0.2), ylim=c(0, 0.2))
abline(0, 1, lty=2)
abline(0, mean(vec_se_nls/vec_se_optim), lty=3)
text(0.025, 0.005, label=mean(vec_se_nls/vec_se_optim))

residuals <- CostFunction(vec_p_optim, vec_t, vec_data, ss=FALSE)
residuals <- sort(residuals)
res_Q <- (seq(length(residuals)) - 0.5)/length(residuals)
print(res_Q)
norm_Q <- qnorm(res_Q, mean=mean(residuals), sd=sd(residuals))

print(residuals)
plot(residuals, norm_Q, xlim=mean(residuals) - c(-3*sd(residuals), 3*sd(residuals)), 
     ylim=mean(residuals) - c(-3*sd(residuals), 3*sd(residuals)))

# segments(0, -3*sd(residuals), 0, 3*sd(residuals), lty=2)

segments(mean(residuals), -3*sd(residuals), mean(residuals), 3*sd(residuals), lty=3)
abline(0, 1)
break

# Plot the data against the underlying model:
# Plot the recovered parameters against the true values:
plot(vec_p, vec_p_optim, xlim=c(0, 5), ylim=c(0, 5))
abline(0, 1, lty=2)


