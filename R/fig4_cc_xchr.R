# fig 4: CC genotype probabilities on X chr

library(qtl2)
library(here)
library(broman)

cachedir <- here("CCapp/_cache")

cc <- readRDS(file.path(cachedir, "cc.rds"))
pmap <- readRDS(file.path(cachedir, "pmap.rds"))
pr_cc <- readRDS(file.path(cachedir, "probs_cc.rds"))
pr_genril <- readRDS(file.path(cachedir, "probs_genril.rds"))
pr_cc038_wrong <- readRDS(file.path(cachedir, "probs_cc038_wrong.rds"))


pdf(here("Figs/fig4_cc_xchr.pdf"), height=7, width=7, pointsize=14)

par(mfrow=c(3,1), mar=c(4.1, 4.1, 2.1, 1.1), las=1)
str <- "CC038/GeniUnc"

par(col.main="darkslateblue")
plot_genoprob(pr_cc, pmap, str, "X", main="more-exact model", yaxt="n",
              mgp.x=c(1.6, 0.2, 0), xlab="Chr X position (Mbp)", xaxt="n")
axis(side=2, at=8:1, names(CCcolors), tick=FALSE, mgp=c(0, 0.3, 0))
xat <- seq(20, 160, by=20)
axis(side=1, xat, mgp=c(1.6, 0.2, 0), tick=FALSE)
u <- par("usr")
text(u[1]-diff(u[1:2])*0.08, u[4]+diff(u[3:4])*0.05, "a", font=2, xpd=TRUE, cex=1.3)

plot_genoprob(pr_genril, pmap, str, "X", main="approximate model", yaxt="n",
              mgp.x=c(1.6, 0.2, 0), xlab="Chr X position (Mbp)", xaxt="n")
axis(side=2, at=8:1, names(CCcolors), tick=FALSE, mgp=c(0, 0.3, 0))
axis(side=1, xat, mgp=c(1.6, 0.2, 0), tick=FALSE)
text(u[1]-diff(u[1:2])*0.08, u[4]+diff(u[3:4])*0.05, "b", font=2, xpd=TRUE, cex=1.3)

plot_genoprob(pr_cc038_wrong, pmap, chr="X", main="wrong cross information", yaxt="n",
              mgp.x=c(1.6, 0.2, 0), xlab="Chr X position (Mbp)", xaxt="n")
axis(side=2, at=8:1, names(CCcolors), tick=FALSE, mgp=c(0, 0.3, 0))
axis(side=1, xat, mgp=c(1.6, 0.2, 0), tick=FALSE)
text(u[1]-diff(u[1:2])*0.08, u[4]+diff(u[3:4])*0.05, "c", font=2, xpd=TRUE, cex=1.3)
dev.off()
