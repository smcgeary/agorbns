
library(deSolve)

graphics.off()


# kon is 3.9e8 M-1 s-1

kons <- c(9, 8, 7, 6, 5, 4, 3, 2, 1, 0.5, 0.01, 0.005, 0.001)

cols <- c("red", rep(c("purple", "forestgreen", "blue", "lightblue"), 3))

# dev.new(xpos=20, ypos=20, height=5, width=3)

# plot(c(0), c(0), type="l", col="white", xlim=c(0, 150), ylim=c(0, 0.5))

k_obs <- c(0, 0, 0, 0)

koff = 0.0035



frac_equil <- rep(0, length(kons))


  Dil <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dx <- l*a*kon - koff*x
      da <- koff*x - l*a*kon
      dl <- koff*x - l*a*kon
      list(c(dl, da, dx))
    })
  }

MonomericModel <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dr    <- -2*kon*r*a + koff*ra
    da    <- -2*kon*r*a + koff*ra
    dra   <- 2*kon*r*a - koff*ra + 2*koff*ra_2 - kon*ra*a
    dra_2 <- kon*ra*a - 2*koff*ra_2
    list(c(dr, da, dra, dra_2))
  })
}

parameters <- c(kon=1*1e8/(1e12)*60, koff=0.05)

A <- 0.5

state <- c(r=0, a=0, ra=5*A, ra_2=A)

times <- seq(0, 100, by=1)

out <- ode(y=state, times=times, func=MonomericModel, parms=parameters, maxsteps=10000)

r   <- out[, 2]
ra  <- out[, 4]
ra2 <- out[, 5]

frac_RNA <- (ra + ra2)/(r + ra + ra2)
frac_site  <- (ra + 2*ra2)/(2*(r + ra + ra2))
frac_RNA_peter <- 2*frac_site - frac_site^2

dev.new(xpos=720, ypos=20, height=4, width=4)

plot(out[, 1], frac_site, type="l",
     ylim=c(0, 1), xlim=c(0, max(times)))

lines(out[, 1], frac_RNA, col="red")
lines(out[, 1], frac_RNA_peter, col="purple")
legend(0, 1, legend=c("frac site", "frac RNA", "frac RNA_peter"),
       col=c("black", "red", "purple"), lwd=1)

break





for (i in 1:length(kons)) {

  parameters <- c(kon = kons[i]*1e9/(1e12)*60, koff=koff)

  A <- 0.5

  state <- c(l=0, a=0, x=A)






  times <- seq(0, 2000, by=0.01)
  Exp <- function(times, b1, b0, k) {
    b1*(1 - exp(-1*k*times)) + b0
  }

  Cost <- function(pars, y) {
    b1 <- pars[1]
    b0 <- pars[2]
    k  <- pars[3]
    y_exp <- Exp(times, b1, b0, k)
    sum((y - y_exp)^2)
  }
  out <- ode(y=state, times=times, func=Dil, parms=parameters, maxsteps=10000)

  pre_factor <- 1/(2*kons[i]*sqrt(koff + 4*A*kons[i]))

  C0 <- koff + 4*A*kons[i]
  C2 <- (koff + 2*A*kons[i])*sqrt(C0)
  C3 <- sqrt(koff)*(koff + 4*A*kons[i])
  # C4 <- tanh(1/2*sqrt(koff)*sqrt(koff + 4*)]

  out_alt <- A^2*kons[i]*times - koff*times^2/2 - A*kons[i]*times^2 + kons[i]*times^3/3

  pars_init <- c(0.25, 0.25, 0.035)
  opt <- optim(pars_init, Cost, y=out[, 4])
  print(opt)
  opt_pars <- opt$par
  b1 <- opt_pars[1]
  b0 <- opt_pars[2]
  k <- opt_pars[3]
  k_obs[i] <- k
  head(out)
  lines(times, out[, 4], col=cols[i])
  frac_equil[i] <- out[nrow(out), 4]/out[1, 4]
  print(length(times))
  y_plot <- Exp(times, b1, b0, k)
  print(head(y_plot))
  lines(times, Exp(times, b1, b0, k), col=cols[i], lty=2)
  lines(times, out_alt, col=cols[i], lty=3, lwd=2)

}

dev.new(height=5, width=5, xpos=320, ypos=20)


plot(k_obs, koff/(1 - frac_equil), col=cols)

koff_calc <- koff/(1 - frac_equil)

lm_fit <- lm(koff_calc ~ k_obs)

coefs <- lm_fit$coefficients
b <- coefs[1]
m <- coefs[2]

segments(x0=0.1, y0=0.1, x1=koff, y1=koff, lty=2)

x_fit <- seq(0, 0.1, length.out=20)
y_fit <- x_fit*m + b

lines(x_fit, y_fit, col="gray")
