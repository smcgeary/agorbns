source("general/general.R")
graphics.off()

poll_data <- data.frame(yes=c(3, 12), no=c(8, 1),
                        row.names=c("idea club", "journal club"))
# poll_data_df <- as.data.frame(poll_data)


z <- 1.96
poll_data$n <- rowSums(poll_data)
poll_data$p <- poll_data$yes/poll_data$n
poll_data$ci <- z*sqrt(poll_data$p*(1 - poll_data$p)/poll_data$n)

xmin <- 0
xmax <- 5
ymin <- 0
ymax <- 1
dev.new(xpos=20, ypos=20, height=5, width=5)
par(kPlotParameters)
BlankPlot()
AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="Yes (%)", percent=TRUE)

text(x=c(0.5, 1.5), 0, labels=c("Idea club", "Journal club"), srt=30,
                     xpd=NA, adj=c(1, 1))
# barplot(poll_data$p)
segments(0, 0.5, 3, 0.5, lty=2)
text(3.2, 0.5, labels="Threshold for official\ninvitation according to\nthe poll", adj=0, xpd=NA)
rect(c(0, 1), c(0, 0), c(1, 2), poll_data$p, col=c("forestgreen", "violet"), border=NA)

arrows(c(0.5, 1.5), poll_data$p - poll_data$ci, c(0.5, 1.5), poll_data$p + poll_data$ci,
     length=0.05*par()$cex, angle=90,
     code=3, xpd=NA)


print(poll_data)
dev.copy2pdf(file="temp_poll.pdf", useDingbats=FALSE)