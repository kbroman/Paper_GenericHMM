# QTL results for Al-Barghouthi et al (2021) data
# show LOD curves + difference
# admit differences from the original paper

library(qtl2)
library(here)

cachedir <- here("DOapp/_cache")
file1 <- file.path(cachedir, "scan1_do.rds")
file2 <- file.path(cachedir, "scan1_genail.rds")
out_scan1_do <- readRDS(file1)
out_scan1_genail <- readRDS(file2)
pmap <- readRDS(file.path(cachedir, "pmap.rds"))

pdf(here("Figs/fig3_do_qtl.pdf"), height=7.5, width=9, pointsize=10)
layout(rbind(1,2), height=c(2,1))
par(mar=c(3.1, 4.1, 2.1, 0.6))

plot(out_scan1_do, pmap, lod="TMD", mgp.x=c(1.8, 0.3, 0), mgp.y=c(2.1, 0.3, 0), lwd=1)
plot(out_scan1_genail, pmap, lod="TMD", add=TRUE, col="violetred", lty=3, lwd=1)
u <- par("usr")
text(u[1]-diff(u[1:2])*0.06, u[4]+diff(u[3:4])*0.03, "a", font=2, cex=1.5, xpd=TRUE)
legend("topright", c("DO", "genAIL8"), lwd=2, lty=c(1,3), col=c("slateblue", "violetred"), bg="gray90")


plot(out_scan1_do - out_scan1_genail, pmap, ylim=c(-0.02, 0.02), lod="TMD", lwd=1,
     ylab="LOD(DO) - LOD(genAIL8)", mgp.x=c(1.8, 0.3, 0), mgp.y=c(2.8, 0.3, 0))
u <- par("usr")
text(u[1]-diff(u[1:2])*0.06, u[4]+diff(u[3:4])*0.14, "b", font=2, cex=1.5, xpd=TRUE)

dev.off()
