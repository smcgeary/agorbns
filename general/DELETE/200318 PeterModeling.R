
library(deSolve)

graphics.off()


x_data <- c(   0,    1,    3,   10,   30,   60,  150)
y_data <- c(0.89, 0.82, 0.77, 0.59, 0.42, 0.48, 0.47)

MonomericModel <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dr    <- -2*kon*r*a + koff*ra
    da    <- -2*kon*r*a + koff*ra - kon*ra*a + 2*koff*ra2
    dra   <-  2*kon*r*a - koff*ra - kon*ra*a + 2*koff*ra2
    dra2  <-                        kon*ra*a - 2*koff*ra2
    list(c(dr, da, dra, dra2))
  })
}


ModelSimulator <- function(times, kon, koff, theta) {
	parameters <- c(kon, koff)
	# Define the total Ago and ligands as 1 pM
	A <- 1
	R <- 0.5
	# Generate the initial conditions
	r_0   <- R*(1 - theta)^2
	ra_0  <- R*2*theta*(1 - theta)
	ra2_0 <- R*theta^2
	a_0 <- A - ra_0 - 2*ra2_0
	state <- c(r=r_0, a=a_0, ra=ra_0, ra2=ra2_0)
	names(state) <- c("r", "a", "ra", "ra2")

	times <- seq(0, max(times), by=1)
	print(times)

	out_model <- ode(y=state, times=times, func=MonomericModel,
	                 parms=parameters, maxsteps=10000)
	out_model
}

CostFunction <- function(pars, x_data, y_data) {
	print("Cost")
	out <- ModelSimulator(times=x_data, kon=10^pars[1], koff=10^pars[2],
	                        theta=1/(1 + exp(-pars[3])))
	out <- out[which(out[, 1] %in% x_data), ]
	t <- out[, 1]
	r <- out[, 2]
	a <- out[, 3]
	ra <- out[, 4]
	ra2 <- out[, 5]
	r_bound <- ra + ra2
	r_tot <- r + ra + ra2
	y_sim <- r_bound/r_tot
	# plot(x_data, y_data, type="p", ylim=c(0, 1))
	# lines(x_data, y_sim, col="blue")
	return(sum((y_sim - y_data)^2))
}

pars_init <- c(kon=log10(1*1e8/(1e12)*60), koff=log10(0.05), theta=0)

opt <- optim(pars_init, CostFunction, x_data=x_data, y_data=y_data)

pars      <- opt$par
kon_fit   <- 10^pars[1]
koff_fit  <- 10^pars[2]
theta_fit <- 1/(1 + exp(-pars[3]))
out <- ModelSimulator(times=x_data, kon=kon_fit, koff=koff_fit,
	                        theta=theta_fit)

t   <- out[, 1]
r   <- out[, 2]
a   <- out[, 3]
ra  <- out[, 4]
ra2 <- out[, 5]

frac_RNA <- (ra + ra2)/(r + ra + ra2)


dev.new(height=4, width=4)
plot(x_data, y_data, type="p", ylim=c(0, 1), xlim=c(0, max(t)))
lines(t, frac_RNA, col="blue")

dev.new(xpos=400, height=4, width=4)
plot(t, a + ra + 2*ra2, ylim=c(0, 1))
lines(t, r + ra + ra2, col="red")