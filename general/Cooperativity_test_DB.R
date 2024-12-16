
# A custom function I made that lets you call functions within functions and not
# have to re-type in the arguments every time.
SubfunctionCall <- function(f, ...) {
  # Check for the supplied arguments in "...":
  f.args.supplied <- list(...)
  # Remove these arguments from the list of all the arguments of function "f":
  f.arg.names <- setdiff(names(formals(f)), names(f.args.supplied))
  # Check for these arguments in the enviromnent:
  f.args <- lapply(f.arg.names, function(arg) {
    arg.temp <- try(dynGet(arg), silent=TRUE)
    if (class(arg.temp) != "try-error" & class(arg.temp) != "error") arg.temp
    else "remove"
  })
  names(f.args) <- f.arg.names
  # Combine the found arguments with the supplied and perform the function:
  f.args.list <- as.list(c(f.args[f.args != "remove"], f.args.supplied))
  do.call(f, f.args.list)
}




# Function that gives the total bound sites, as a function of 
# K1 : Kd of AGO-miR-1
# K2 : Kd of AGO-miR-124
# delta : The ratio describing the fold-increase in the Kd of double-site
# 	binding in comparison to single-site binding.
# m1 : The concentration of free AGO-miR-1
# m2 : The concentration of free AGO-miR-124
FractionalOccupancy <- function(K1, K2, delta, m1, m2) {
	(K2*m1 + K1*m2 + 2*delta*m2*m1)/(2*(K1*K2 + K2*m1 + K1*m2 + delta*m1*m2))
}


# Function that gives total amount of AGO-miR-1 and AGO-miR-124, as a function
#	of
# K1 : Kd of AGO-miR-1
# K2 : Kd of AGO-miR-124
# delta : The ratio describing the fold-increase in the Kd of double-site
# 	binding in comparison to single-site binding.
# m1 : The concentration of free AGO-miR-1
# m2 : The concentration of free AGO-miR-124
TotalmiRNAs <- function(K1, K2, delta, m1, m2) {
	x0 <- K1*K2/(K1*K2 + K2*m1 + K1*m2 + delta*m1*m2)
	x1 <- m1*x0/K1
	x2 <- m2*x0/K2
	x3 <- m1*m2*x0/(K1*K2/delta)
	m1_T <- m1 + x1 + x3
	m2_T <- m2 + x2 + x3
	return(c(m1_T, m2_T))
}

# Function that quantifies the degree of error between the known total AGO-miR-1
# AGO-miR-124 concentrations, and the values given by the current estimates of
# the free AGO-miR-1 and AGO-miR-124. This is a function of
# K1 : Kd of AGO-miR-1
# K2 : Kd of AGO-miR-124
# delta : The ratio describing the fold-increase in the Kd of double-site
# 	binding in comparison to single-site binding.
# m1_T : The known concentration of total AGO-miR-1
# m2_T : The known concentration of free AGO-miR-124
ResidualFunction <- function(par, K1, K2, delta, m1_T, m2_T) {
	# Pulling the two free miRNA concentrations out of the 'par' variable.
	m1 <- par[1]
	m2 <- par[2]
	# Calculate the total AGO-miR-1 and AGO-miR-124 as a function of K1, K2,
	# delta, m1, and m2.
	m_T <- SubfunctionCall(TotalmiRNAs)
	# Return the summed squared residuals between the known and estimated total
	# AGO-miR-1 and AGO-miR-124.
	 sum((m_T - c(m1_T, m2_T))^2)
}


# Function that returns the concentrations of free AGO-miR-1 and AGO-miR-124,
# as a function of
# K1 : Kd of AGO-miR-1
# K2 : Kd of AGO-miR-124
# delta : The ratio describing the fold-increase in the Kd of double-site
# 	binding in comparison to single-site binding.
# m1_T : The known concentration of total AGO-miR-1
# m2_T : The known concentration of free AGO-miR-124
GetFreemiRNAs <- function(K1, K2, delta, m1_T, m2_T) {
	# Initialize the starting guess for the concentrations of free AGO-miR-1 and
	# AGO-miR-124, at 50% of their total concentrations.
	pars <- c(m1_T, m2_T)/2
	# optimization routine that starts with the paramter values 'pars', and
	# arrives at the correct values by minimizing the output of the above error
	# function.
	opt <- optim(par=pars, fn=ResidualFunction, K1=K1, K2=K2, delta=delta,
				 m1_T=m1_T, m2_T=m2_T, method="L-BFGS-B",
	             lower=c(0, 0), upper=c(m1_T, m2_T))
	print(opt)
	return(opt$par)
}



# These remaining functions directly calculate the values of interest as a
# function of the variables that we in principle know (K1, K2, delta, m1_T, and
# m2_T)

# This calculates the total occupancy across both sites.
FullModel <- function(K1, K2, delta, m1_T, m2_T) {
	m12 <- SubfunctionCall(GetFreemiRNAs)
	m1 <- m12[1]
	m2 <- m12[2]
	return(SubfunctionCall(FractionalOccupancy))
}

# This calculates the concentration of the unbound target RNA.
X0 <- function(K1, K2, delta, m1_T, m2_T) {
	m12 <- SubfunctionCall(GetFreemiRNAs)
	m1 <- m12[1]
	m2 <- m12[2]
	K1*K2/(K1*K2 + K2*m1 + K1*m2 + delta*m1*m2)
}

# This calculates the concentration of singly bound target RNA (by either 
# AGO-miR-1 or AGO-miR-124).
SingleBound <- function(K1, K2, delta, m1_T, m2_T) {
	m12 <- SubfunctionCall(GetFreemiRNAs)
	m1 <- m12[1]
	m2 <- m12[2]
	x0 <- K1*K2/(K1*K2 + K2*m1 + K1*m2 + delta*m1*m2)
	m1*x0/K1 + m2*x0/K2
}

# This calculates the concentration of the doubly bound target RNA.
DoubleBound <- function(K1, K2, delta, m1_T, m2_T) {
	m12 <- SubfunctionCall(GetFreemiRNAs)
	m1 <- m12[1]
	m2 <- m12[2]
	x0 <- K1*K2/(K1*K2 + K2*m1 + K1*m2 + delta*m1*m2)
	m1*m2*x0/(K1*K2/delta)
}


# Concentration series approximating the experiments on your sheet.
m1_T_series <- c(0.3, 1.5, 3, 6, 15, 30, 60, 150, 300, 600)
m2_T <- 30

# Chosen for visual agreement with the figures on your sheet.
K1 <- 30
K2 <- 50
delta <- 1

# Get the binding occupancy across the AGO-miR-1 concentration series:
occs <- sapply(m1_T_series, FullModel, K1=K1, K2=K2, delta=delta, m2_T=m2_T)
# Get the free target RNA concentrations across the AGO-miR-1 concentration 
# series:
x0 <- sapply(m1_T_series, X0, K1=K1, K2=K2, delta=delta, m2_T=m2_T)
# Get the singly bound target RNA concentration across the AGO-miR-1
# concentration series:
x1.2 <- sapply(m1_T_series, SingleBound, K1=K1, K2=K2, delta=delta, m2_T=m2_T)
# Get the doubly bound target RNA concentration across the AGO-miR-1
# concentration series:
x3 <- sapply(m1_T_series, DoubleBound, K1=K1, K2=K2, delta=delta, m2_T=m2_T)

dev.new(xpos=20, ypos=20, height=5, width=5)
plot(m1_T_series, occs, log='x')

dev.new(xpos=520, ypos=20, height=5, width=5)

plot(m1_T_series, x1.2, ylim=c(0, 1), log='x', col="red")
legend("topright", col=c("red", "blue", "green"), pch=20, legend=c("single", "double", "free"))
points(m1_T_series, x3, col="blue")
points(m1_T_series, x0, col="green")
points(m1_T_series, (x0 + x1.2 + x3), col="black")