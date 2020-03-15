mdsPlot<-function (dat, numPositions = 1000, sampNames = NULL, sampGroups = NULL, 
          xlim, ylim, pch = 16, pal = brewer.pal(8, "Dark2"), cex=1.25,
          legendPos = "bottomleft", legendNCol=1, main = NULL){
  if (is(dat, "MethylSet") || is(dat, "RGChannelSet")) {
    b <- getBeta(dat)
  }
  else if (is(dat, "matrix")) {
    b <- dat
  }
  else {
    stop("dat must be an 'MethylSet', 'RGChannelSet', or 'matrix'.")
  }
  if (is.null(main)) {
    main <- sprintf("Beta MDS\n%d most variable positions", numPositions)
  }
  o <- order(rowVars(b), decreasing = TRUE)[seq_len(numPositions)]
  d <- dist(t(b[o, ]))
  fit <- cmdscale(d)
  if (missing(xlim)) 
    xlim <- range(fit[, 1]) * 1.2
  if (missing(ylim)) 
    ylim <- range(fit[, 2]) * 1.2
  if (is.null(sampGroups)) 
  sampGroups <- rep(1, numPositions)
  sampGroups <- as.factor(sampGroups)
  col <- pal[sampGroups]
  if (is.null(sampNames)) {
    plot(x = fit[, 1], y = fit[, 2], col = col, pch = pch, xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = main,cex=cex)
  }
  else {
    plot(x = 0, y = 0, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = main,cex=cex)
    text(x = fit[, 1], y = fit[, 2], sampNames, col = col)
  }
  numGroups <- length(levels(sampGroups))
  if (missing(legendNCol)) 
    legendNCol <- numGroups
  if (numGroups > 1) {
    legend(x = legendPos, legend = levels(sampGroups), ncol = legendNCol, text.col = pal[seq_len(numGroups)],bty="n")
  }
}
