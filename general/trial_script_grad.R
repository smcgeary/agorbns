source("general/general.R")

sXc <- SitesXCounts("miR-1", sitelist="paper")


InitializeEquilSitePars <- function(sXc, combined=TRUE) {
  l <- SubfunctionCall(GetInputEquil)
  # print(l)
  data <- SubfunctionCall(GetDataEquil)
  # print(data)
  kds <- log10(Norm(l)/Norm(rowSums(data, -ncol(data)))/
               max(Norm(l)/Norm(rowSums(data[,-ncol(data)]))))
  # kds <- kds*0
  kds <- kds - kds[length(kds)]
  names(kds) <- paste0(rownames(sXc), "_Kd")
  pars <- c(kds, bg=-4, AGO=1)
  return(pars)
}


pars <- InitializeEquilSitePars(sXc)


print(CostNew(pars, sXc))
print(CostC(pars, sXc))
plot(GradientNew(pars, sXc), GradientC(pars, sXc))
print(cbind(GradientNew(pars, sXc), GradientC(pars, sXc)))
abline(0, 1)
break

# l <- GetInputEquil(sXc)
# kds <-10^pars[1:length(l)]
# A <- 10^pars[length(l)+2]
# bg <- 10^pars[length(l)+1]

# time1 <- proc.time()[3]
# a_R <- FreeAgo(A=A, kds=kds, l=l)
# time2 <- proc.time()[3]
# a_C <- FreeAgoC(A=A, kds=kds, l=l)
# time3 <- proc.time()[3]
# print(a_R)
# print(time2 - time1)
# print(a_C)
# print(time3 - time2)
# print(a_R/a_C)

time1 <- proc.time()[3]
out <- GradientC(pars, sXc)
time2 <- proc.time()[3]
print("C done")
print(time2 - time1)
time2 <- proc.time()[3]
out.n <- GradientNew(pars, sXc)
time3 <- proc.time()[3]
print("R done")
print(time3 - time2)

out.num <- grad(CostNew, pars, sXc=sXc)
print(cbind(out, out.n, out.num))

print(CostNew(pars, sXc))
print(CostC(pars, sXc))


break
