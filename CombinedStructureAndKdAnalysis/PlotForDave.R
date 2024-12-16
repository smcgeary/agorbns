source("general/general.R")

ind_8mer <- which(I.data$site=="8mer")
ind_6mer <- which(I.data$site == "6mer" & I.data$n > 0)
dev.new(xpos = 20, ypos = 20)


    xmin <- 10
    xmax <- 200
    ymin <- 10
    ymax <- 200


plot(c.final[ind_8mer,column]*sum(data[,column])/I.data$n[ind_8mer],
	data[ind_8mer,column]/I.data$n[ind_8mer],
	log='xy',
	xlim = c(xmin, xmax),
	ylim = c(ymin, ymax),
	col = sapply(as.character(I.data$flank[ind_8mer]),GetColorFunction),
	lwd=2,axes=FALSE,ann=FALSE)
   



    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))

    xs <- unique(c(xmin,xs[xs >= xmin & xs <= xmax],xmax))
    ys <- unique(c(ymin,ys[ys >= ymin & ys <= ymax],ymax))

    # xmin <- min(xs)
    # xmax <- max(xs)
    # ymin <- min(ys)
    # ymax <- max(ys)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
    yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))


   axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=ymin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=yl,
         labels=sapply(yl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=ys, labels=FALSE,
         pos=xmin, lwd=2)

        title("8mer", font.main = 1, cex.main = 1.5, line=-2, adj=0.1)
        cor_text <- round(
                      cor(log(data[ind_8mer,column]/I.data$n[ind_8mer]),
	log(c.final[ind_8mer,column]*sum(data[,column])/I.data$n[ind_8mer])),
                      digits = 3
                    )

        title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)
        title(xlab = "Predicted enrichment from\ncomputational folding", cex.lab = 1.5)
        title(ylab = "Experimental enrichment", cex.lab =1.5, line = 1.5)
dev.new(xpos = 1000, ypos = 20)

    xmin <- 1
    xmax <- 40
    ymin <- 1
    ymax <- 40



plot(c.final[ind_6mer,column]*sum(data[,column])/I.data$n[ind_6mer],
	data[ind_6mer,column]/I.data$n[ind_6mer],
	log='xy',
	xlim = c(xmin, xmax),
	ylim = c(ymin, ymax),
	col = sapply(as.character(I.data$flank[ind_6mer]),GetColorFunction),
	lwd=2,
	axes=FALSE,
	ann=FALSE)




    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))

    xs <- unique(c(xmin,xs[xs >= xmin & xs <= xmax],xmax))
    ys <- unique(c(ymin,ys[ys >= ymin & ys <= ymax],ymax))

    # xmin <- min(xs)
    # xmax <- max(xs)
    # ymin <- min(ys)
    # ymax <- max(ys)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
    yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))




   axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=ymin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=yl,
         labels=sapply(yl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=ys, labels=FALSE,
         pos=xmin, lwd=2)

        title("6mer", font.main = 1, cex.main = 1.5, line=-2, adj=0.1)
        cor_text <- round(
                      cor(log(data[ind_6mer,column]/I.data$n[ind_6mer]),
	log(c.final[ind_6mer,column]*sum(data[,column])/I.data$n[ind_6mer])),
                      digits = 3
                    )

        title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)
                title(xlab = "Predicted enrichment from\ncomputational folding", cex.lab = 1.5)
        title(ylab = "Experimental enrichment", cex.lab =1.5, line = 1.5)
