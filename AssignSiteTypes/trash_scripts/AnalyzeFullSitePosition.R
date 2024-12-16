################################################################################
#GenerateSiteTypeKds.py
################################################################################

# # Initial parameters and constants.

site_cols <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_colors.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}



args = commandArgs(trailingOnly=TRUE)
mirna = args[1]
start = as.numeric(args[2])
stop  = as.numeric(args[3])

conditions <- c("I", "0", "0.4", "1.26", "4", "12.6", "40")
cond_cols <- c("black", "gray", "red", "orange", "green", "blue", "purple")

stockago <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/SolveForKds/k_c_stockago.txt",
  row.names=1, header=FALSE, sep="\t")

k.c.stockago = stockago[mirna,1]

# Get the names of all the site types.
sites <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.",
                     mirna, ".txt"),stringsAsFactors=FALSE)[,1]
print(sites)
names(cond_cols) <- conditions
MakeEnsembleFigures <- function(mirna, condition) {
  sites.pos5p.name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/equilibrium/site_position5p_counts/", condition, "_",
                   start, "-", stop, ".txt")
  site.pos.5p <- read.table(sites.pos5p.name, row.names=1, header=TRUE)
  sites.pos3p.name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/equilibrium/site_position3p_counts/", condition, "_",
                   start, "-", stop, ".txt")
  site.pos.3p <- read.table(sites.pos3p.name, row.names=1, header=TRUE)
  centered_sites <- c("11mer-m3.13", "12mer-m3.14",
                      "11mer-m4.14", "12mer-m4.15")

  centered_index <- sapply(centered_sites, function(site){
                             which(rownames(site.pos.5p) == site)
                           }
                           )

  site.pos.5p <- site.pos.5p[-centered_index,]
  site.pos.3p <- site.pos.3p[-centered_index,]


  # centered.pos.5p <- colSums(site.pos.5p[centered_index,])

  # centered.pos.3p <- colSums(site.pos.3p[centered_index,])

  # site.pos.5p <- rbind(site.pos.5p[-centered_index,],centered.pos.5p)
  # rownames(site.pos.5p)[length(rownames(site.pos.5p))] <- "Centered"
  # site.pos.3p <- rbind(site.pos.3p[-centered_index,],centered.pos.3p)
  # rownames(site.pos.3p)[length(rownames(site.pos.3p))] <- "Centered"
  setEPS(width = 10)

  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/figures/site_positions/", mirna,"/by_condition/",condition,"_",
                   start, "-", stop, ".eps"))
  par(par_plots)

  plot(x    = 1: dim(site.pos.5p)[2],
       y    = rep(0, dim(site.pos.5p)[2]),
       xlim = c(1, 2*dim(site.pos.5p)[2]),
       ylim = c(0, 1),
       col  = "white")
    
  sapply(rownames(site.pos.5p), function(name) {
            temp5 <- site.pos.5p[name,]
            temp3 <- site.pos.3p[name,]
            diff <- c(sapply(1:length(temp5), function(ind){
            if (temp5[ind] != 0) {
            return(c(unlist(which(temp3[1,]==c(temp5[ind]))))-ind)

            }
            }))
            offset <- getmode(unlist(diff))[[1]]
            x = seq(1,(37+start+stop-offset))+offset/2
            y = Cumul(site.pos.5p[name,1:(dim(site.pos.5p)[2]-offset)])
           lines(x = seq(1,(37+start+stop-offset))+offset/2,
                 y = Cumul(site.pos.5p[name,1:(dim(site.pos.5p)[2]-offset)]),
                 col = site_cols[name, 1])
    })

  legend("topright",
         legend = rownames(site.pos.3p),
         lwd = 2,
         col = sapply(rownames(site.pos.3p), function(name){
          site_cols[name,1]
          }))
  dev.off()
}
sapply(conditions, function(condition){
  MakeEnsembleFigures(mirna,condition)
  } )

MakeSiteFigures <- function(mirna, site) {
  outs <- t(sapply(conditions, function(condition) {
      sites.pos5p.name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/equilibrium/site_position5p_counts/", condition, "_",
                   start, "-", stop, ".txt")
      return(read.table(sites.pos5p.name, row.names=1, header=TRUE)[site,])
  }))
  outs.plot <- outs
  setEPS()

  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/figures/site_positions/", mirna,"/by_site/",site,"_",
                   start, "-", stop, ".eps"))
  par(par_plots)
  x = seq(dim(outs.plot)[2])
  plot(x    = x,
       y    = rep(0, dim(outs.plot)[2]),
       xlim = c(0, max(x)+1),
       ylim = c(0, 1),
       col  = "white",
       axes = FALSE)
  title(main = mirna, line=0.5, adj=0.1)
  print(site_cols)
  print(site)
  print(site_cols[site,])
  title(main = site, col.main=site_cols[site,], line=-1, adj=0.1)

  title(xlab = "Site position")
  title(ylab = "Cumulative Distribution")
  axis(1, at=seq(0, max(x)+1, by=max(1, floor(max(x) / 20))), pos=0, lwd=2,
       labels=FALSE, tck=-0.01)
  axis(1, at=seq(0, max(x)+1, by=max(1, floor(max(x) / 20))*5), pos=0, lwd=2)
  axis(2, at=c(0,0.2,0.4,0.6,0.8,1),
       pos=0, lwd=2, hadj=0.8)

  legend(x=1,
         y=0.98,
         legend = conditions,
         lwd = 2,
         cex = 1.2,
         col = cond_cols,
         bty = "n")

  sapply(rownames(outs.plot), function(name) {
           lines(x = c(0,rep(1:dim(outs.plot)[2],each=2),dim(outs.plot)[2]+1),
                 y = rep(c(0,Cumul(unlist(outs.plot[name,]))),each=2),
                 col = cond_cols[name])
    })
  dev.off()

}
sapply(sites, function(site) {
  print(site)
  MakeSiteFigures(mirna,site)
})


