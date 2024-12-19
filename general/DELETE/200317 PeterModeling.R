
library(deSolve)

graphics.off()

MonomericModel <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dr    <- -2*kon*r*a + koff*ra
    da    <- -2*kon*r*a + koff*ra
    dra   <- 2*kon*r*a - koff*ra + 2*koff*ra2 - kon*ra*a
    dra2 <- kon*ra*a - 2*koff*ra2
    list(c(dr, da, dra, dra2))
  })
}

MonomericModelGood <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dr    <- -2*kon*r*a + koff*ra
    da    <- -(2*r + ra)*kon*a + koff*(ra + 2*ra2)
    dra   <- 2*kon*r*a - koff*ra + 2*koff*ra2 - kon*ra*a
    dra2 <- kon*ra*a - 2*koff*ra2
    list(c(dr, da, dra, dra2))
  })
}


parameters <- c(kon=1*1e8/(1e12)*60, koff=0.05)

ra2_initial <- 0.5

state <- c(r=0, a=0, ra=0, ra2=ra2_initial)

times <- seq(0, 100, by=1)

out <- ode(y=state, times=times, func=MonomericModel, parms=parameters, maxsteps=10000)
out2 <- ode(y=state, times=times, func=MonomericModelGood, parms=parameters, maxsteps=10000)

t   <- out[, 1]
r   <- out[, 2]
a   <- out[, 3]
ra  <- out[, 4]
ra2 <- out[, 5]

A_tot <- a + ra + 2*ra2

R_tot <- r + ra + ra2

frac_site  <- (ra + 2*ra2)/(2*(r + ra + ra2))
frac_RNA <- (ra + ra2)/(r + ra + ra2)
frac_RNA_peter <- 2*frac_site - frac_site^2

dev.new(height=4, width=4)
plot(out[, 1], frac_site, type="l", ylim=c(0, 1), xlim=c(0, max(times)))
lines(out[, 1], frac_RNA, col="red")
lines(out[, 1], frac_RNA_peter, col="purple")

dev.new(xpos=400, height=4, width=4)

plot(t, A_tot, type="l")
lines(t, R_tot, col="red")


out <- ode(y=state, times=times, func=MonomericModelGood, parms=parameters, maxsteps=10000)

t   <- out[, 1]
r   <- out[, 2]
a   <- out[, 3]
ra  <- out[, 4]
ra2 <- out[, 5]

A_tot <- a + ra + 2*ra2

R_tot <- r + ra + ra2

frac_site  <- (ra + 2*ra2)/(2*(r + ra + ra2))
frac_RNA <- (ra + ra2)/(r + ra + ra2)
frac_RNA_peter <- 2*frac_site - frac_site^2

# dev.new(height=4, width=4)
dev.set(2)
lines(out[, 1], frac_site, lty=2)
lines(out[, 1], frac_RNA, col="red", lty=2)
lines(out[, 1], frac_RNA_peter, col="purple", lty=3, lwd=2)

dev.set(3)

lines(t, A_tot, lwd=2)
lines(t, R_tot, col="red", lwd=3, lty=2)


break

legend(0, 1, legend=c("frac site", "frac RNA", "frac RNA_peter"),
       col=c("black", "red", "purple"), lwd=1)

dev.new(xpos=400, height=4, width=4)

plot(t, A_tot, type="l")
lines(t, R_tot, col="red")

