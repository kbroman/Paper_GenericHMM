# illustration of genome reconstruction in MPP
# use DO data from Al-Bargouthi et al. (2021)
# mouse #42, chr 12, 30-33 Mbp

library(qtl2)
library(qtl2fst)
library(broman)
library(here)

hom <- cumsum(1:8)
pr <- readRDS(here("DOapp/_cache/probs_do/probs_do_fstindex.rds"))
load(here("DOapp/_cache/cross_basic_cleaned.Rdata"))
pmap <- readRDS(here("DOapp/_cache/pmap.rds"))

d <- 1
chr <- 14
pos <- 101.143898832531
ind <- "49"
map <- pmap[[chr]]
interval <- pos+(d/2)*c(-1,1)
map <- map[map >= interval[1] & map <= interval[2]]
mar<- names(map)
nmar <- length(mar)

p <- pr[[14]][ind,,mar]

# adjust map so no small gaps
dmap <- diff(map)
dmap[dmap < 0.01] <- 0.01
amap <- setNames(cumsum(c(map[1], dmap)), names(map))

fg <- cross_basic$founder_geno[[chr]][,mar]
g <- cross_basic$geno[[chr]][ind,mar]
fg[fg==0] <- NA
g[g==0] <- NA

# swap alleles so 1 is more common
afreq <- colMeans(fg==1, na.rm=TRUE)
fg[,afreq < 0.5] <- 4-fg[,afreq<0.5]
g[afreq < 0.5] <- 4-g[afreq<0.5]

# swap D <-> F so I can use blue/yellow and still conform to official CC colors
fg[c(1,4,6),] <- fg[c(6,1,4),]
rn <- rownames(p)
rn_new <- gsub("x", "F", gsub("F", "A", gsub("A", "D", gsub("D", "x", rn))))
rn_new <- sapply(lapply(strsplit(rn_new, ""), sort), paste, collapse="")
rownames(p) <- rn_new
p <- p[rn,]


nf <- nrow(fg)

pdf(here("Figs/fig1_genome_reconstr.pdf"), height=6.5, width=6, pointsize=14)

layout(rbind(1,2,3,4), height=c(6,1,3.5,2.5))
par(mar=c(0.6, 4.1, 1.6, 0.3))

linecolor <- "gray90"
pointcolor <- c("white", "gray", "black")
ticks <- seq(100.6, 101.6, by=0.2)
probcolors <- brocolors("web")[c("blue", "green", "yellow")]
probs <- c("DD", "AD")


# A: founder genotypes
plot(0,0,type="n", xlab="", ylab="Founders", xaxt="n", yaxt="n",
     xlim=range(amap), ylim=c(nf + 0.5, 0.5), yaxs="i", mgp=c(2.1, 0, 0))
abline(h=1:nf, col=linecolor)
axis(side=2, at=1:nf, labels=LETTERS[1:nf], tick=FALSE, mgp=c(0, 0.3, 0), las=1)
abline(v=ticks, col=linecolor)
for(i in 1:nf) {
    points(amap, rep(i, nmar), pch=21, bg=pointcolor[fg[i,]])
}
box()
u <- par("usr")
text(u[1]-diff(u[1:2])*0.10, u[4]+diff(u[3:4])*0.04, "a", font=2, xpd=TRUE, cex=1.5)


# B: genotypes of DO individual
plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n",
     xlim=range(amap), ylim=c(0.5, 1.5), yaxs="i", mgp=c(2.1, 0, 0))
axis(side=2, at=1, label=paste0("DO-", ind), tick=FALSE, las=1, mgp=c(0, 0.2, 0))
abline(v=ticks, col=linecolor)
abline(h=1:nf, col=linecolor)
points(amap, rep(1, nmar), pch=21, bg=pointcolor[g])
box()
u <- par("usr")
text(u[1]-diff(u[1:2])*0.10, u[4]+diff(u[3:4])*0.50, "b", font=2, xpd=TRUE, cex=1.5)


# C: genotype probabilities
plot(0,0,type="n", xlab="", ylab="genotype probability", xaxt="n", yaxt="n",
     xlim=range(amap), ylim=c(-0.02, 1.02), yaxs="i", mgp=c(2.1, 0.3, 0))
axis(side=2, at=seq(0, 1, by=0.2), mgp=c(0, 0.2, 0), las=1, tick=FALSE)
abline(h=seq(0, 1, by=0.2), col=linecolor)
abline(v=ticks, col=linecolor)
amap2 <- amap; amap2[1] <- amap[1]-1; amap2[length(amap)] <- amap2[length(amap)]+1
u <- par("usr")
for(i in seq_along(probs)) {
    lines(amap2, p[probs[i],], col=probcolors[i], lwd=2)
    text(u[1]+diff(u[1:2])*c(0.02, 0.98)[i], 0.9, probs[i], col=probcolors[i])
}
box()
text(u[1]-diff(u[1:2])*0.10, u[4]+diff(u[3:4])*0.06, "c", font=2, xpd=TRUE, cex=1.5)


# D: inferred genotypes
par(mar=c(3.1, 4.1, 1.6, 0.3))
plot(0,0,type="n", xlab=paste("Chr", chr, "position (Mbp)"), ylab="", yaxt="n", xaxt="n",
     xlim=range(amap), ylim=c(0.5, 2), yaxs="i", mgp=c(1.8, 0.6, 0))
abline(v=ticks, col=linecolor)
axis(side=1, at=ticks, tick=FALSE, mgp=c(0, 0.3, 0))
axis(side=2, at=1.25, label=paste0("DO-", ind), tick=FALSE, las=1, mgp=c(0, 0.2, 0))
u <- par("usr")
text(u[1]-diff(u[1:2])*0.10, u[4]+diff(u[3:4])*0.12, "d", font=2, xpd=TRUE, cex=1.5)

rect(u[1], 0.7, u[2], 1.1, border=probcolors[1], col=probcolors[1])
text(u[1]+diff(u[1:2])*0.02, 0.9, col="white", "D")

left <- max(amap2[p[probs[1],] > 0.9])
right <- min(amap2[p[probs[2],] > 0.9])
rect(u[1], 1.3, left, 1.7, border=probcolors[1], col=probcolors[1])
text(u[1]+diff(u[1:2])*0.02, 1.5, col="white", "D")
rect(right, 1.3, u[2], 1.7, border=probcolors[3], col=probcolors[3])
text(u[1]+diff(u[1:2])*0.98, 1.5, col="black", "A")

box()

dev.off()
