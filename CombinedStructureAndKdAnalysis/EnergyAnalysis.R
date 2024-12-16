
site_cols <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_info.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 1, drop = FALSE]

graphics.off()
GetKdTable <- function(mirna,experiment) {
	stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")

	k.c.stockago <- stockago[mirna, experiment]
	table <- read.table(file = paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
		       mirna, "/",experiment,"/kds_with_structure/",
		       "-3--3_1-15_preliminary_dinucs_protein.txt"))
	table <- t(table)
	rownames(table) <- table[, 1]
	table <- table[, -1]
	kd_length <- (nrow(table) - 1) 
	return(table[1:kd_length,])
}

kds_sm <- GetKdTable(mirna, "equilibrium")
kds_nb <- GetKdTable(mirna, "equilibrium_nb")
kds_nb_seed <- GetKdTable(mirna, "equil_seed_nb")
cols = c(site_cols[rownames(kds_nb),][1:243],rep("chartreuse",15),"darkslateblue",rep("gray",5))

dev.new(xpos = 20, ypos = 20, height = 10, width = 10)
par(mfrow = c(3,3))
plot(kds_sm[, 1],
	 kds_nb[, 1],
	 col = cols,
	 lwd = 2,
	 pch = 19,
	 xlim = c(-6, 2),
	 ylim = c(-6, 2))
plot(kds_sm[,ncol(kds_sm)],
	 kds_nb[,ncol(kds_nb)],
	 col = cols,
	 lwd = 2,
	 pch = 19,
	 xlim = c(-6, 2),
	 ylim = c(-6, 2))
plot(seq(1,min(ncol(kds_nb),ncol(kds_sm)), length = 50),
	 sapply(seq(1,min(ncol(kds_nb),ncol(kds_sm)), length = 50),
	 	    function(i){cor(as.numeric(kds_sm[,ncol(kds_sm)-min(ncol(kds_nb),ncol(kds_sm)) + i]),as.numeric(kds_nb[,ncol(kds_sm)-min(ncol(kds_nb),ncol(kds_sm)) + i]))}),
	 ann = FALSE)

plot(kds_sm[, 1],
	 kds_nb_seed[, 1],
	 col = cols,
	 lwd = 2,
	 pch = 19,
	 xlim = c(-6, 2),
	 ylim = c(-6, 2))
plot(kds_sm[,ncol(kds_sm)],
	 kds_nb_seed[,ncol(kds_nb_seed)],
	 col = cols,
	 lwd = 2,
	 pch = 19,
	 xlim = c(-6, 2),
	 ylim = c(-6, 2))
plot(seq(1,min(ncol(kds_nb_seed),ncol(kds_sm)), length = 50),
	 sapply(seq(1,min(ncol(kds_nb_seed),ncol(kds_sm)), length = 50),
	 	    function(i){cor(as.numeric(kds_sm[,i]),as.numeric(kds_nb_seed[,i]))}),
	 ann = FALSE)

plot(kds_nb[, 1],
	 kds_nb_seed[, 1],
	 col = cols,
	 lwd = 2,
	 pch = 19,
	 xlim = c(-6, 2),
	 ylim = c(-6, 2))
plot(kds_nb[,ncol(kds_nb)],
	 kds_nb_seed[,ncol(kds_nb_seed)],
	 col = cols,
	 lwd = 2,
	 pch = 19,
	 xlim = c(-6, 2),
	 ylim = c(-6, 2))
plot(seq(1,min(ncol(kds_nb_seed),ncol(kds_nb)), length = 50),
	 sapply(seq(1,min(ncol(kds_nb_seed),ncol(kds_nb)), length = 50),
	 	    function(i){cor(as.numeric(kds_nb_seed[,i]),as.numeric(kds_nb[,i]))}),
	 ann = FALSE)




dinuc_freqs <- cbind(as.numeric(kds_sm[c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG"),
	                                   ncol(kds_sm)]),
	                 as.numeric(kds_nb[c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG"),
	                                   ncol(kds_nb)]))
names(dinuc_freqs) <- c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG")



all_6mers <- c("6mer-A1",
			   "6mer",
			   "6mer-m8",
			   "6mer-m4.9",
			   "6mer-m5.10",
			   "6mer-m6.11",
			   "6mer-m7.12",
			   "6mer-m8.13",
			   "6mer-m9.14",
			   "6mer-m10.15",
			   "6mer-m11.16",
			   "6mer-m12.17",
			   "6mer-m13.18",
			   "6mer-m14.19",
			   "6mer-m15.20",
			   "6mer-m16.21",
			   "6mer-m17.22")

all_7mers <- c("7mer-A1",
			   "7mer-m8",
			   "7mer-m3.9",
			   "7mer-m4.10",
			   "7mer-m5.11",
			   "7mer-m6.12",
			   "7mer-m7.13",
			   "7mer-m8.14",
			   "7mer-m9.15",
			   "7mer-m10.16",
			   "7mer-m11.17",
			   "7mer-m12.18",
			   "7mer-m13.19",
			   "7mer-m14.20",
			   "7mer-m15.21",
			   "7mer-m16.22")

all_8mers <- c("8mer",
			   "8mer-m2.9",
			   "8mer-m3.10",
			   "8mer-m4.11",
			   "8mer-m5.12",
			   "8mer-m6.13",
			   "8mer-m7.14",
			   "8mer-m8.15",
			   "8mer-m9.16",
			   "8mer-m10.17",
			   "8mer-m11.18",
			   "8mer-m12.19",
			   "8mer-m13.20",
			   "8mer-m14.21",
			   "8mer-m15.22")

all_9mers <- c("9mer-m1.9",
			   "9mer-m2.10",
			   "9mer-m3.11",
			   "9mer-m4.12",
			   "9mer-m5.13",
			   "9mer-m6.14",
			   "9mer-m7.15",
			   "9mer-m8.16",
			   "9mer-m9.17",
			   "9mer-m10.18",
			   "9mer-m11.19",
			   "9mer-m12.20",
			   "9mer-m13.21",
			   "9mer-m14.22")


all_10mers <- c("10mer-m1.10",
			    "10mer-m2.11",
			    "10mer-m3.12",
			    "10mer-m4.13",
			    "10mer-m5.14",
			    "10mer-m6.15",
			    "10mer-m7.16",
			    "10mer-m8.17",
			    "10mer-m9.18",
			    "10mer-m10.19",
			    "10mer-m11.20",
			    "10mer-m12.21",
			    "10mer-m13.22")

all_11mers <- c("11mer-m1.11",
			    "11mer-m2.12",
			    "11mer-m5.15",
			    "11mer-m6.16",
			    "11mer-m7.17",
			    "11mer-m8.18",
			    "11mer-m9.19",
			    "11mer-m10.20",
			    "11mer-m11.21",
			    "11mer-m12.22")

all_12mers <- c("12mer-m1.12",
			    "12mer-m2.13",
			    "12mer-m5.16",
			    "12mer-m6.17",
			    "12mer-m7.18",
			    "12mer-m8.19",
			    "12mer-m9.20",
			    "12mer-m10.21",
			    "12mer-m11.22")

all_13mers <- c("13mer-m1.13",
			    "13mer-m2.14",
			    "13mer-m4.16",
			    "13mer-m5.17",
			    "13mer-m6.18",
			    "13mer-m7.19",
			    "13mer-m8.20",
			    "13mer-m10.22")

kds_sm_final <- kds_sm[,ncol(kds_sm), drop = FALSE]
kds_nb_final <- kds_nb[, ncol(kds_nb),drop = FALSE]
kds_nb_seed_final <- kds_nb_seed[, ncol(kds_nb_seed),drop = FALSE]

kd_6mer_profile <- cbind(kds_sm_final[all_6mers, ],
						 kds_nb_final[all_6mers, ],
						 kds_nb_seed_final[all_6mers, ])

kd_7mer_profile <- cbind(kds_sm_final[all_7mers, ],
						 kds_nb_final[all_7mers, ],
						 kds_nb_seed_final[all_7mers, ])

kd_8mer_profile <- cbind(kds_sm_final[all_8mers, ],
						 kds_nb_final[all_8mers, ],
						 kds_nb_seed_final[all_8mers, ])

kd_9mer_profile <- cbind(kds_sm_final[all_9mers, ],
						 kds_nb_final[all_9mers, ],
						 kds_nb_seed_final[all_9mers, ])

kd_10mer_profile <- cbind(kds_sm_final[all_10mers, ],
						 kds_nb_final[all_10mers, ],
						 kds_nb_seed_final[all_10mers, ])

kd_11mer_profile <- cbind(kds_sm_final[all_11mers, ],
						 kds_nb_final[all_11mers, ],
						 kds_nb_seed_final[all_11mers, ])

kd_12mer_profile <- cbind(kds_sm_final[all_12mers, ],
						 kds_nb_final[all_12mers, ],
						 kds_nb_seed_final[all_12mers, ])

kd_13mer_profile <- cbind(kds_sm_final[all_13mers, ],
						 kds_nb_final[all_13mers, ],
						 kds_nb_seed_final[all_13mers, ])



# kd_7mer_profile <- sapply(c("equilibrium", "equilibrium_nb"),
# 					   function(experiment) {
# 					   	kds <- GetKdTable("let-7a", experiment)
# 					   	return(kds[all_7mers, ncol(kds)])
# 					   })

# kd_8mer_profile <- sapply(c("equilibrium", "equilibrium_nb"),
# 					   function(experiment) {
# 					   	kds <- GetKdTable("let-7a", experiment)
# 					   	return(kds[all_8mers, ncol(kds)])
# 					   })
# kd_9mer_profile <- sapply(c("equilibrium", "equilibrium_nb"),
# 					   function(experiment) {
# 					   	kds <- GetKdTable("let-7a", experiment)
# 					   	return(kds[all_9mers, ncol(kds)])
# 					   })
# kd_10mer_profile <- sapply(c("equilibrium", "equilibrium_nb"),
# 					   function(experiment) {
# 					   	kds <- GetKdTable("let-7a", experiment)
# 					   	return(kds[all_10mers, ncol(kds)])
# 					   })

# kd_11mer_profile <- sapply(c("equilibrium", "equilibrium_nb"),
# 					   function(experiment) {
# 					   	kds <- GetKdTable("let-7a", experiment)
# 					   	return(kds[all_11mers, ncol(kds)])
# 					   })

# kd_12mer_profile <- sapply(c("equilibrium", "equilibrium_nb"),
# 					   function(experiment) {
# 					   	kds <- GetKdTable("let-7a", experiment)
# 					   	return(kds[all_12mers, ncol(kds)])
# 					   })

# kd_13mer_profile <- sapply(c("equilibrium", "equilibrium_nb"),
# 					   function(experiment) {
# 					   	kds <- GetKdTable("let-7a", experiment)
# 					   	return(kds[all_13mers, ncol(kds)])
# 					   })

colors = c("forestgreen", "blue", "magenta", "red", "orange")
dev.new(xpos = 20, ypos = 20, width = 16, height = 10)
par(mfrow=c(2,4))
plot(seq(1:length(all_6mers)),kd_6mer_profile[,2],type="o",col="white",ylim=c(-5,3))
sapply(1:3,function(col) {
	lines(seq(1:length(all_6mers)),
		  kd_6mer_profile[,col],
		  type="o",col=colors[col],
		  lwd = 2)
	})
legend("topleft",c("SM","NB", "NB_seed"),
	   lwd = 2,
	   col = colors)

plot(seq(1:length(all_7mers)),kd_7mer_profile[,2],type="o",col="white",ylim=c(-5,3))
sapply(1:3,function(col) {
	lines(seq(1:length(all_7mers)),
		  kd_7mer_profile[,col],
		  type="o",col=colors[col],
		  lwd = 2)
	})
legend("topleft",c("SM","NB", "NB_seed"),
	   lwd = 2,
	   col = colors)

plot(seq(1:length(all_8mers)),kd_8mer_profile[,2],type="o",col="white",ylim=c(-5,3))
sapply(1:3,function(col) {
	lines(seq(1:length(all_8mers)),
		  kd_8mer_profile[,col],
		  type="o",col=colors[col],
		  lwd = 2)
	})
legend("topleft",c("SM","NB", "NB_seed"),
	   lwd = 2,
	   col = colors)

plot(seq(1:length(all_9mers)),kd_9mer_profile[,2],type="o",col="white",ylim=c(-5,3))
sapply(1:3,function(col) {
	lines(seq(1:length(all_9mers)),
		  kd_9mer_profile[,col],
		  type="o",col=colors[col],
		  lwd = 2)
	})
legend("topleft",c("SM","NB", "NB_seed"),
	   lwd = 2,
	   col = colors)
plot(seq(1:length(all_10mers)),kd_10mer_profile[,2],type="o",col="white",ylim=c(-6,3))
sapply(1:3,function(col) {
	lines(seq(1:length(all_10mers)),
		  kd_10mer_profile[,col],
		  type="o",col=colors[col],
		  lwd = 2)
	})
legend("topleft",c("SM","NB", "NB_seed"),
	   lwd = 2,
	   col = colors)

plot(seq(1,length(all_11mers)),kd_11mer_profile[,2],type="o",col="white",ylim=c(-6,3))
sapply(1:3,function(col) {
	lines(seq(1,length(all_11mers)),
		  kd_11mer_profile[,col],
		  type="o",col=colors[col],
		  lwd = 2)
	})
legend("topleft",c("SM","NB", "NB_seed"),
	   lwd = 2,
	   col = colors)

plot(seq(1,length(all_12mers)),kd_12mer_profile[,2],type="o",col="white",ylim=c(-6,3))
sapply(1:3,function(col) {
	lines(seq(1,length(all_12mers)),
		  kd_12mer_profile[,col],
		  type="o",col=colors[col],
		  lwd = 2)
	})
legend("topleft",c("SM","NB", "NB_seed"),
	   lwd = 2,
	   col = colors)

plot(seq(1,length(all_13mers)),kd_13mer_profile[,2],type="o",col="white",ylim=c(-6,3))
sapply(1:3,function(col) {
	lines(seq(1,length(all_13mers)),
		  kd_13mer_profile[,col],
		  type="o",col=colors[col],
		  lwd = 2)
	})
legend("topleft",c("SM","NB", "NB_seed"),
	   lwd = 2,
	   col = colors)







