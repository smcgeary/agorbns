# site <- "GCACTTTA"

whole_thing <- function(site){
nb_pos_I    <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_nb/site_positions/",site,"/I_5-5_current.txt"))
nb_pos_0    <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_nb/site_positions/",site,"/0_5-5_current.txt"))
nb_pos_0.4  <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_nb/site_positions/",site,"/0.4_5-5_current.txt"))
nb_pos_1.26 <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_nb/site_positions/",site,"/1.26_5-5_current.txt"))
nb_pos_4    <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_nb/site_positions/",site,"/4_5-5_current.txt"))
nb_pos_12.6 <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_nb/site_positions/",site,"/12.6_5-5_current.txt"))
nb_pos_40   <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_nb/site_positions/",site,"/40_5-5_current.txt"))


sm_pos_I    <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/site_positions/",site,"/I_5-5_current.txt"))
sm_pos_0    <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/site_positions/",site,"/0_5-5_current.txt"))
sm_pos_0.4  <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/site_positions/",site,"/0.4_5-5_current.txt"))
sm_pos_1.26 <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/site_positions/",site,"/1.26_5-5_current.txt"))
sm_pos_4    <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/site_positions/",site,"/4_5-5_current.txt"))
sm_pos_12.6 <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/site_positions/",site,"/12.6_5-5_current.txt"))
sm_pos_40   <- read.table(file=paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/site_positions/",site,"/40_5-5_current.txt"))


par(mfrow = c(1, 2))



x_positions <- rep(rownames(nb_pos_I), each=2)
x_positions <- x_positions[-length(x_positions)]
y_func <- function(data) {
	out <- sapply(seq(1,nrow(nb_pos_I)),function(i) {sum(Norm(data[,1])[1:(i+1)])})
	out <- rep(out, each =2)
	return(out[-1])
}

plot(x_positions, y_func(nb_pos_I), type = "l", ylim = c(0, 1), lwd = 2)
lines(x_positions,y_func(nb_pos_0),    lwd = 2, col = "gray")
lines(x_positions,y_func(nb_pos_0.4),  lwd = 2, col = "red")
lines(x_positions,y_func(nb_pos_1.26), lwd = 2, col = "orangered")
lines(x_positions,y_func(nb_pos_4),    lwd = 2, col = "green")
lines(x_positions,y_func(nb_pos_12.6), lwd = 2, col = "blue")
lines(x_positions,y_func(nb_pos_40),   lwd = 2, col = "purple")
plot(x_positions, y_func(sm_pos_I), type = "l", ylim = c(0, 1), lwd = 2)
lines(x_positions,y_func(sm_pos_0),    lwd = 2, col = "gray")
lines(x_positions,y_func(sm_pos_0.4),  lwd = 2, col = "red")
lines(x_positions,y_func(sm_pos_1.26), lwd = 2, col = "orangered")
lines(x_positions,y_func(sm_pos_4),    lwd = 2, col = "green")
lines(x_positions,y_func(sm_pos_12.6), lwd = 2, col = "blue")
lines(x_positions,y_func(sm_pos_40),   lwd = 2, col = "purple")
}

