# plot of transition probability on autosome and X chromosome
# with "do" cross type, "genail8" cross type, and the difference
# use like r = 0.01 or 0.05 and generation 8,16,24,32

library(qtl2)
library(broman)
library(here)

r <- seq(0, 0.5, len=501)
gen <- c(8, 16, 24, 32)

# X chromosome (male)
doX <- lapply(gen, function(g) qtl2:::test_stepmatrix("do", r, TRUE, FALSE, g))
genailX <- lapply(gen, function(g) qtl2:::test_stepmatrix("genail8", r, TRUE, FALSE, c(g+5, rep(1, 8))))

# autosome
doA <- lapply(gen, function(g) qtl2:::test_stepmatrix("do", r, FALSE, FALSE, g))
genailA <- lapply(gen, function(g) qtl2:::test_stepmatrix("genail8", r, FALSE, FALSE, c(g+5, rep(1, 8))))

# pull out just the first AA->AA case
trprob2matrix <-
function(v)
{
    ngen <- length(v)
    nrec <- length(v[[1]])
    result <- matrix(nrow=nrec, ncol=ngen)

    for(g in 1:ngen) {
        for(r in 1:nrec) {
                result[r,g] <- v[[g]][[r]][1,1]
        }
    }

    exp(result)
}

doX <- 1-trprob2matrix(doX)
genailX <- 1-trprob2matrix(genailX)

doA <- 1-sqrt(trprob2matrix(doA))
genailA <- 1-sqrt(trprob2matrix(genailA))

linecolors <- brocolors("web")[c("blue", "red", "green", "orange")]

pdf(here("Figs/fig2_do_trmatrix.pdf"), height=5, width=6.5, pointsize=10)
par(mfrow=c(2,2), mar=c(3.1, 3.1, 2.1, 1.1))

grayplot(r, doA[,1], type="n", ylim=c(0, 0.88), yaxs="i", xaxs="i",
         xlab="recombination fraction", ylab="exchange probability",
         main="Autosome", mgp.x=c(1.6, 0.3, 0), mgp.y=c(1.8, 0.2, 0))
for(i in 1:ncol(doA)) {
    lines(r, doA[,i], col=linecolors[i])
    lines(r, genailA[,i], col=linecolors[i], lty=2)
}
u <- par("usr")
text(u[1]-diff(u[1:2])*0.12, u[4]+diff(u[3:4])*0.08, "a", font=2, cex=1.5, xpd=TRUE)

grayplot(r, doA[,1]/genailA[,1], type="n", ylim=c(0.99, 1.05), yaxs="i", xaxs="i",
         xlab="recombination fraction", ylab="ratio of exchange probabilities",
         main="Autosome", mgp.x=c(1.6, 0.3, 0), mgp.y=c(2.2, 0.2, 0))
for(i in 1:ncol(doA)) {
    lines(r, doA[,i]/genailA[,i], col=linecolors[i])
}
u <- par("usr")
text(u[1]-diff(u[1:2])*0.12, u[4]+diff(u[3:4])*0.08, "c", font=2, cex=1.5, xpd=TRUE)

legend("topright", paste("gen", gen), col=linecolors, lwd=2, bg="gray90")



grayplot(r, doX[,1], type="n", ylim=c(0, 0.88), yaxs="i", xaxs="i",
         xlab="recombination fraction", ylab="exchange probability",
         main="X chromosome", mgp.x=c(1.6, 0.3, 0), mgp.y=c(1.8, 0.2, 0))
for(i in 1:ncol(doX)) {
    lines(r, doX[,i], col=linecolors[i])
    lines(r, genailX[,i], col=linecolors[i], lty=2)
}
u <- par("usr")
text(u[1]-diff(u[1:2])*0.12, u[4]+diff(u[3:4])*0.08, "b", font=2, cex=1.5, xpd=TRUE)

grayplot(r, doX[,1]/genailX[,1], type="n", ylim=c(0.99, 1.05), yaxs="i", xaxs="i",
         xlab="recombination fraction", ylab="ratio of exchange probabilities",
         main="X chromosome", mgp.x=c(1.6, 0.3, 0), mgp.y=c(2.2, 0.2, 0))
for(i in 1:ncol(doX)) {
    lines(r, doX[,i]/genailX[,i], col=linecolors[i])
}
u <- par("usr")
text(u[1]-diff(u[1:2])*0.12, u[4]+diff(u[3:4])*0.08, "d", font=2, cex=1.5, xpd=TRUE)

dev.off()
