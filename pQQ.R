# this is for qqplot
pQQ<-function (pvals, nlabs = 6, conf = 0.95, lim, mark = 0.05, ...) 
{
    if (any(is.na(pvals))) 
        stop("Missing p-values not allowed", call. = F)
    sim <- F
    .pvals <- sort(pvals)
    .logpvals <- -log10(.pvals)
    .order <- seq(along = .pvals)
    .n <- length(.pvals)
    .aux <- .order/(.n + 1)
    .logaux <- -log10(.aux)
    .xt <- 1/2
    .fix.order <- c(.xt, .order)
    .fix.logaux <- -log10(.fix.order/(.n + 1))
    .b.median <- qbeta(p = 1/2, shape1 = .fix.order, shape2 = .n - 
        .fix.order + 1)
    .f.median <- approxfun(x = .fix.logaux, y = -log10(.b.median))
    .adj.logaux <- .f.median(.logaux)
    .adj.fix.logaux <- .f.median(.fix.logaux)
    if (missing(lim)) {
        .lim <- c(0, max(.adj.logaux, .logpvals)) * 1.05
    }
    else .lim <- lim
    plot(.logaux, .logaux, type = "n", pch = 16, cex = 0.7, xlim = .lim, 
        ylim = .lim, xlab = "Expected P-value (-log10 scale)", 
        ylab = "Observed P-value (-log10 scale)", font = 2, lwd = 2, 
        font.lab = 2, bty = "l", xaxs = "i", yaxs = "i", axes = T, 
        las = 1, pty = "s", ...)
    if (is.numeric(conf)) {
        .b.lower <- qbeta(p = (1 - conf)/2, shape1 = .fix.order, 
            shape2 = .n - .fix.order + 1)
        .b.upper <- qbeta(p = 1 - (1 - conf)/2, shape1 = .fix.order, 
            shape2 = .n - .fix.order + 1)
        .sm.lower <- approx(x = .adj.fix.logaux, y = -log10(.b.lower), 
            xout = seq(.lim[1], .lim[2], length.out = 1000))
        .sm.upper <- approx(x = .adj.fix.logaux, y = -log10(.b.upper), 
            xout = seq(.lim[1], .lim[2], length.out = 1000))
        segments(x0 = .sm.lower$x, y0 = .sm.lower$y, x1 = .sm.upper$x, 
            y1 = .sm.upper$y, col = "grey")
        lines(.adj.fix.logaux, -log10(.b.lower), col = "black", 
            lwd = 1)
        lines(.adj.fix.logaux, -log10(.b.upper), col = "black", 
            lwd = 1)
        if (sim) {
            .conf <- f.QQconf(nGenes = length(.pvals), quantiles = c((1 - 
                conf)/2, 1 - (1 - conf)/2))
            .s.lower <- rev(.conf$quant[1, ])
            .s.upper <- rev(.conf$quant[2, ])
            lines(.adj.logaux, .s.lower, lty = 2, col = "red", 
                lwd = 2)
            lines(.adj.logaux, .s.upper, lty = 2, col = "red", 
                lwd = 2)
        }
    }
    points(.adj.logaux, .logpvals, pch = 16, cex = 0.7)
    box(lwd = 2, bty = "l")
    abline(a = 0, b = 1, lwd = 2, col = "red")
    if (is.numeric(mark)) {
        psig <- -log10(mark)
        lines(c(psig, psig), c(0, psig), lty = 3, lwd = 1.4)
        lines(c(0, psig), c(psig, psig), lty = 3, lwd = 1.4)
    }
    if (nlabs > 0) {
        .ind <- seq(length.out = nlabs)
        text(.adj.logaux[.ind], .logpvals[.ind] - 0.05, names(.pvals)[.ind], 
            srt = 270, font = 2, cex = 0.7, adj = 0, col = "black", 
            xaxs = "i", yaxs = "i")
    }
    return(invisible())
}
