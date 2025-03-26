plotfgsea <- function(pathway, stats,fgseaRes, gene.set.name, class.name, posClass, negClass,
                      gene.set.title.name = 'Gene Set'){
  stopifnot(!missing(pathway))
  stopifnot(!missing(stats))
  stopifnot(!missing(fgseaRes))
  stopifnot(gene.set.name %in% fgseaRes$pathway )
  stats = sort(stats, decreasing = T)
  metric.range <- c(min(stats), max(stats))
  gsea.enrichment.score= round(fgseaRes$ES[which(fgseaRes$pathway==gene.set.name)],2)
  gsea.normalized.enrichment.score= round(fgseaRes$NES[which(fgseaRes$pathway==gene.set.name)],2)
  gsea.p.value = fgseaRes$pval[which(fgseaRes$pathway==gene.set.name)]
  gsea.p.value = formatC(gsea.p.value, format = "e", digits = 2)
  gsea.fdr = fgseaRes$padj[which(fgseaRes$pathway==gene.set.name)]
  gsea.fdr = formatC(gsea.fdr, format = "e", digits = 2)

  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * abs(statsAdj)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  enrichment.score.range <- c(min(toPlot$y), max(toPlot$y))

  
  
  
  gsea.layout <- layout(matrix(c(1, 2, 3)), heights = c(1.7, 0.3, 0.2))

  par(mar = c(0, 3, 6, 3))
  
  
  plot(toPlot$x, toPlot$y, type = "l",
       col = "#A4C361", lwd = 2, lty=1, pch=19, 
       cex = 2, cex.axis = 2, 
       xaxt = "n",xaxs = "i", xlab = "", ylab = "Enrichment score",
       ylim = enrichment.score.range,
       main = list(gene.set.title.name, 
                   font = 1, cex = 2),
       
       panel.first = {
         abline(h = seq(round(enrichment.score.range[1], digits = 1),
                        enrichment.score.range[2], 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot.coordinates <- par("usr")
  if(gsea.enrichment.score < 0) {
    text(length(stats) * 0.01, plot.coordinates[3] * 0.98, cex = 2,
         paste( "P:", gsea.p.value,#"FDR:", gsea.fdr,#
                "\nNES:",gsea.normalized.enrichment.score, "\n"), adj = c(0, 0))
  } else {
    text(length(stats) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03), cex = 2,
         paste( "P:", gsea.p.value,#"FDR:", gsea.fdr,#
                "\nNES:",gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
  }

  par(mar = c(0, 3, 0, 3))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(stats)))
  abline(v = toPlot$x, lwd = 0.75)

  par(mar = c(1, 3, 0, 3))
  rank.colors <- stats - metric.range[1]
  rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
  rank.colors <- ceiling(rank.colors * 255 + 1)
  tryCatch({
    rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
  }, error = function(e) {
    stop("Please use the metric.range argument to provide a metric range that",
         "includes all metric values")
  })
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(stats)))
  box()
  
  
  text(length(stats) / 2, 0.7, cex = 2,
       labels = ifelse(!missing(class.name), class.name, ''))
  text(length(stats) * 0.01, 0.7, cex = 2,
       ifelse(!missing(posClass), posClass, 'Positive'), adj = c(0, NA))
  text(length(stats) * 0.99, 0.7, cex = 2,
       ifelse(!missing(negClass), negClass, 'Negative'), adj = c(1, NA))
}

