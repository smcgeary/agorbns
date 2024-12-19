
library(deSolve)

TEST_SYSTEM <- function(t, state, parameters) {
	with(as.list(c(state, parameters)), {
		dX0 <- a - (B0 + k0) * X0
		dX1 <- k0*X0 - B1*X1
		list(c(dX0, dX1))
	})
}
parameters_dependent <- c(a = 1e3, B0 = 0.001, k0 = 0.06, B1=0.2)
parameters_independent <- c(a = 1.5e3, B0 = 0.070, k0 = 0.021, B1=0.070)

state <- c(X0=0, X1=0)
times <- seq(0, 150, length.out=1000)
out_dependent <- ode(y = state, times = times, func = TEST_SYSTEM, parms = parameters_dependent)
out_independent <- ode(y = state, times = times, func = TEST_SYSTEM, parms = parameters_independent)

plot(out_dependent[, 1], out_dependent[, 2], log='y', type="l", ylim=c(10, 1e5))
lines(out_dependent[, 1], out_dependent[, 3], col="red")

lines(out_independent[, 1], out_independent[, 2], lty=2, lwd=2)
lines(out_independent[, 1], out_independent[, 3], lty=2, lwd=2, col="red")

legend("bottomright", legend=c("miRNA, dep. pars", "tailed, dep. pars",
                               "miRNA, ind. pars", "tailed, ind. pars"),
	   lty=c(1, 1, 2, 2), lwd=c(1, 1, 2, 2),
	   col=c("black", "red", "black", "red"))


