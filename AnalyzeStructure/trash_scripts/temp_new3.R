ecdf_dens <- function(density){x <- density$x; y <- sapply(seq(1,length(density$y)),function(i){sum(density$y[1:i])/length(x)}); return(data.frame(x,y))}
ecdf_fun <- function(y){return(sapply(seq(1,length(y)),function(i){return(sum(y[1:i])/sum(y))}))}
library(truncnorm)
# D <- density(I^(1/9),n=5000,from=0,to=1)
# D_A <- density(A^(1/9),n=20000,from=0.000002,to=1)
x <- seq(0.00002,1,length=4000)
dots <- seq(1,5000,length=20)
n=9
kd = 0.1
Ago=1
get_beta_params <- function(func) {
	u <- mean(func)
	var <- var(func)
	alpha <- max(((1 - u)/ var - 1 / u)*u^2,1)
	beta <- max(alpha * (1 / u - 1),1)
	return(c(alpha,beta))

}
plot_beta <- function(func,col="black"){
	params <- get_beta_params(func)
	alpha <- params[1]
	beta <- params[2]
	lines(x,pbeta(x,alpha,beta),lwd=2,col=col)
}
plot(x[seq(1,4000,length=20)],ecdf(I^(1/9))(seq(0,1,length=4000))[seq(1,4000,length=20)],type="p",lwd=2,col="blue",xlim=c(0,1),ylim=c(0,1))
points(x[seq(1,4000,length=20)],ecdf(A^(1/9))(seq(0,1,length=4000))[seq(1,4000,length=20)],lwd=2,col="red")
# plot_beta(I^(1/9),col="blue")
# plot_beta(A^(1/9),col="red")
# plot_beta(A_1^(1/9),col="orangered")

# plot_beta(A_4^(1/9),col="purple")
# plot_beta(A_12^(1/9),col="deeppink")

# plot_beta(A_40^(1/9),col="forestgreen")

points(x[seq(1,4000,length=20)],ecdf(A_4^(1/9))(seq(0,1,length=4000))[seq(1,4000,length=20)],lwd=2,col="purple")
points(x[seq(1,4000,length=20)],ecdf(A_40^(1/9))(seq(0,1,length=4000))[seq(1,4000,length=20)],lwd=2,col="cyan")
# points(x[seq(1,4000,length=20)],ecdf(A_12^(1/9))(seq(0,1,length=4000))[seq(1,4000,length=20)],lwd=2,col="deeppink")
# points(x[seq(1,4000,length=20)],ecdf(A_1^(1/9))(seq(0,1,length=4000))[seq(1,4000,length=20)],lwd=2,col="orangered")

Get_flank_parameters <- function(flank){
	data <- I[which(flanks_temp==flank)]
	return(get_beta_params(data^(1/9)))
}
fit_flank_parameters <- function(flank){
	data <- I[which(flanks_temp==flank)]
	data_plot <- ecdf(data^(1/9))(x)
	plot_beta(data,col="gray")
	# solution <- optim(c(0,0),function(pars){
	# 	mean = 1/(1+exp(-pars[1]))
	# 	sd = 10^pars[2]
	# 	data <- ecdf((I^(1/9))[which(flanks_temp==flank)])(x)
	# 	model <- ecdf_fun(dtruncnorm(x,a=0,b=1,mean=mean,sd=sd))
	# 	lines(x,model,lwd=0.5)

	# 	return(sum( (data - model)^2))
	# 	})
	# pars <- solution$par
	# mean = 1/(1+exp(-pars[1]))
	# sd = 10^pars[2]
	# model <- ecdf_fun(dtruncnorm(x,a=0,b=1,mean=mean,sd=sd))
	# lines(x,model,col="red",lwd=2)
	# return(solution)

}
plot_flank_parameters <- function(flank,data,col="gray"){
	params <- Get_flank_parameters(flank)
	data <- I[which(flanks_temp==flank)]
	data_plot <- ecdf(data^(1/9))(x)
	plot(x,data_plot,type="l")
	lines(x,pbeta(x,params[1],params[2]),lwd=2,col=col)


}


get_flank_correlation_window <- function(data_I,data_A){
	data <- I[which(flanks_temp==flank)]
	return(get_beta_params(data^(1/9)))

}
sapply()
# break
# lines(x,ecdf_fun(dtruncnorm(x,a=0,b=1,mean=0.59,sd=0.42)),lwd=3,lty=2,col="red")
data_A <- ecdf(A^(1/9))(x)[seq(1,4000,length=20)]
data_A_4 <- ecdf(A_4^(1/9))(x)[seq(1,4000,length=20)]
data_I <- ecdf(I^(1/9))(x)[seq(1,4000,length=20)]
data_A_40 <- ecdf(A_40^(1/9))(x)[seq(1,4000,length=20)]

# lines(x,ecdf_fun(dtruncnorm(x,a=0,b=1,mean=0.59,sd=0.42)*occ(n,kd,Ago)),lwd=3,lty=2,col="blue")

 source("AnalyzeStructure/temp_new5.R")

input_beta_dist <- matrix(apply(parameters,2,function(col){col[3]*dbeta(x,col[1],col[2])}),nrow=dim(parameters)[2],byrow=TRUE)
cor_ <- function(pars){

	n <- pars[1]
	kd <- 10^pars[2]
	Ago <- 10^pars[3]
	output_beta_dist <- matrix(apply(input_beta_dist,1,function(row){row*occ(n,kd,Ago)}),ncol=dim(input_beta_dist)[2],byrow=TRUE)
	plot(log(total_flanks_A),log(rowSums(output_beta_dist)))
	rss <- 1 - cor(log(total_flanks_A), log(rowSums(output_beta_dist)))^2
	print(rss)
	return(rss)

}
R_ <- function(pars){

	print(pars)
	alpha = 10^pars[1]
	print(alpha)
	beta = 10^pars[2]
	print(beta)
	n <- 1.5
	print(n)
	kd <- 10^pars[3]
	print(kd)
	Ago <- 10^pars[4]
	print(Ago)
	# Ago2 <- 10^pars[6]
	# print(Ago2)

	model_I <- ecdf_fun(dbeta(x,alpha,beta))
	model_A <- ecdf_fun(dbeta(x,alpha,beta)*occ(n,kd,Ago))
	model_A_4 <- ecdf_fun(dbeta(x,alpha,beta)*occ(n,kd,Ago*10))
	model_A_40 <- ecdf_fun(dbeta(x,alpha,beta)*occ(n,kd,Ago*100))


	lines(x,model_I,col=rgb(0,0,1,alpha=0.3))
	lines(x,model_A,col=rgb(1,0,0,alpha=0.3))
	lines(x,model_A_4,col=rgb(0.5,0,0.5,alpha=0.3))
	lines(x,model_A_40,col=rgb(0,1,1,alpha=0.3))
	model_I <- model_I[seq(1,4000,length=20)]
	model_A <- model_A[seq(1,4000,length=20)]
	model_A_4 <- model_A_4[seq(1,4000,length=20)]
	model_A_40 <- model_A_40[seq(1,4000,length=20)]


	rss <- sum((data_A - model_A)^2)+sum((data_I - model_I)^2) + sum((data_A_4 - model_A_4)^2) + sum((data_A_40 - model_A_40)^2)
	print(rss)
	return(rss)

}
break
pars_3 <- optim(pars_3$par,R_,control=list(maxit=20000))
# lines(x,ecdf_fun(dtruncnorm(x^(1/n),a=0,b=1,mean=0.59,sd=0.42)),lwd=2,lty=2,col="blue")
# points(x[dots],ecdf_fun(dtruncnorm(x^(1/n),a=0,b=1,mean=0.59,sd=0.42))[dots],lwd=2,lty=2,col="blue")

# lines(x,ecdf_fun(dtruncnorm(x^(1/n),a=0,b=1,mean=0.59,sd=0.42)*((1/n)*x^((1/n)-1))),lwd=2,lty=2,col="purple")
# points(x[dots],ecdf_fun(dtruncnorm(x^(1/n),a=0.0004,b=1,mean=0.59,sd=0.42)*(1/n)*x^((1/n)-1))[dots],lwd=2,lty=2,col="purple")


# lines(x^(1/9),ecdf(dtruncnorm(x,a=0,b=1,mean=0.59,sd=0.42)),lwd=2,lty=2,col="green")
# lines((x^(1/9))[dots],eddf(dtruncnorm(x,a=0,b=1,mean=0.59,sd=0.42))[dots],lwd=2,lty=2,col="green")
break
dots <- seq(1,5000,length=20)
x <- D$x
points(x[dots],ecdf_dens(D)$y[dots],pch=20)
plot(ecdf(I),add=TRUE,col="orange")
# plot(ecdf(I^(1/9)),add=TRUE,col="orange")
lines(ecdf_dens(D_A),col="red")
points(x[dots],ecdf_dens(D_A)$y[dots],col="red",pch=20)
# y <- ecdf_fun(D_I$y*occ(x^9,1,400))
x_new <- x^9
D_new <- sapply(x_new,function(i){D$y[which(abs(x-i)==min(abs(x-i)))]})
x_div <- (1/9)*x^(-8/9)
lines(x,ecdf_fun(D_new),col="purple")
points(x[dots],ecdf_fun(D_new)[dots],col="purple",pch=20)

lines(x,ecdf_fun(D_new*x_div),col="forestgreen")
points(x[dots],ecdf_fun(D_new*x_div)[dots],col="forestgreen",pch=20)


lines(x_new,ecdf_fun((D$y*x_div)[low_end]),col="blue")
points(x_new[dots],ecdf_fun((D$y*x_div)[low_end])[dots],col="blue",pch=20)

# lines(x^(1/2),D$y*0.2,col="blue")
# x_new <- x^(1/2)
# x_div <- abs(2*x_new)
# points((x^(1/2))[seq(1,10000,length=20)],(D$y*0.2)[seq(1,10000,length=20)],col="blue",pch=20)
# lines(x_new,x_div)
# lines(x_new,D$y*x_div,col="purple")

# lines(ecdf_dens(density(A^(1/9),n=1000,from=0,to=1)),col="red")
# plot(ecdf_dens(density(I)))
# x <- density(I^(1/9),n=1000,from=0,to=1)$x