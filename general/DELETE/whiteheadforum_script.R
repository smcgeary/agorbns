
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")

radiography <- read.csv(paste0("../../Presentations/seminars/170508 Whitehead",
                               " Forum/miR-1_autoradiography.csv"),
                        stringsAsFactors=FALSE)

plot(radiography[1 : 5, 1] / 100,
     radiography[1 : 5, 6],
     log='xy',
     axes=FALSE,
     ann=FALSE)
  # title(main=mirna, line=-1, adj=0.1)
  # title(main=site, col.main=kSiteColors[site,], line=-2.5, adj=0.1)

  x <- radiography[1:6,1]/100*stockago["miR-1","equilibrium"]
  y <- radiography[1:6,6]
  sdev <- radiography[1:6,7]
  sites.norm <- Norm(c.I.tots)
  data.norm <- t(t(data)/colSums(data))


  xmin <- 0
  xmax <- 1
  ymin <- 0
  ymax <- .01
  xs <- seq(xmin, xmax, length = 21)
  xl <- seq(xmin, xmax, length = 5)
  ys <- seq(ymin, ymax, length = 21)
  yl <- seq(ymin, ymax, length = 5)

  plot(x, y,xlim=c(xmin, xmax), ylim=c(ymin, ymax), type="p",
  	pch=19, cex=1.2,
     axes=FALSE, ann=FALSE)        
  # Generate tickmarks for axis.

  axis(1, at=xl,
       labels=sprintf("%.2f",round(xl,2)),
       pos=ymin, lwd=0, cex.axis=1.4, padj=0.2)
  axis(1, at=xs, labels=FALSE,
       pos=ymin, lwd=2)
  # Label the axis at each order of magnitude.

  axis(2, at=yl,
       labels=sprintf("%.2f",round(yl*100,2)),
       pos=xmin, las=2, lwd=0, cex.axis=1.4,hadj=0.9)
  axis(2, at=ys, labels=FALSE,
       pos=xmin, lwd=2)
  fit_model <- lm(y~x)
  r_2 <- round(summary(fit_model)$r.squared, 3)
  title(main = bquote(italic(r)^2 ==  .(r_2)), font.main=1, cex.main=1.4, line=-2, adj=0.9)
  title(xlab = "[AGO2-miR-1] (nM)", cex.lab=1.5, line=2, adj=0.5)
  title(ylab = "% Library", cex.lab=1.4, line=2.5)

	lines(x,lm(y~x)$fitted.values)
dev.copy2pdf(file = "../../Presentations/seminars/170508 Whitehead Forum/170509 Linear_model_radioactivity.pdf")