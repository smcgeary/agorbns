x <- seq(0,100,length=10)
ExpoModel <- function(B,A,k,noise) {
	return(exp(A - k*x) + B + rnorm(length(x),mean=0,sd=noise))
}

par(mfrow=c(1,2))
plot(x,ExpoModel(.1,0.5,0.05,0.0001),ylim=c(0,3),type="l")
data <- ExpoModel(.1,0.5,0.05,0.05)

points(x,data)

FitModelLog <- function(data,pars_init=c(10,0,0)){
	pars_out <- optim(pars_init,function(pars) {
		B <- min(data)/(1 + exp(-pars[1]))
		A <- 10^pars[2]
		k <- 10^pars[3]
		print(pars)
		print(c(B,A,k))
		data_norm <- log(data - B)
		model <- A - k*x
		lines(x,exp(model) + B,col=rgb(0,0,1,alpha=0.1))
		rss <- sum((data_norm - model)^2)
		print(rss)
		return(rss)	
	},control=list(maxit=80000))$par
	pars_out[1] <- min(data)/(1 + exp(-pars_out[1]))
	return(pars_out)
}

FitModelLin <- function(data,pars_init=c(10,0,0)){
	pars_out <- optim(pars_init,function(pars) {
		B <- min(data)/(1 + exp(-pars[1]))
		A <- 10^pars[2]
		k <- 10^pars[3]
		print(pars)
		print(c(B,A,k))
		model <- exp(A - k*x) + B
		lines(x,model,col=rgb(1,0,0,alpha=0.1))
		rss <- sum((data - model)^2)
		print(rss)
		return(rss)	
	},control=list(maxit=80000))$par
	pars_out[1] <- min(data)/(1 + exp(-pars_out[1]))
	return(pars_out)
}
pars_fit_log <- FitModelLog(data)
pars_fit_lin <- FitModelLin(data)
B_fit_log <- pars_fit_log[1]
A_fit_log <- 10^pars_fit_log[2]
k_fit_log <- 10^pars_fit_log[3]

B_fit_lin <- pars_fit_lin[1]
A_fit_lin <- 10^pars_fit_lin[2]
k_fit_lin <- 10^pars_fit_lin[3]



plot(x,ExpoModel(.1,0.5,0.05,0.0001),ylim=c(0,3),type="l",col="black")
points(x,data)
lines(x,exp(A_fit_log - k_fit_log*x)+B_fit_log,col="blue")
lines(x,exp(A_fit_lin - k_fit_lin*x)+B_fit_lin,col="red")

