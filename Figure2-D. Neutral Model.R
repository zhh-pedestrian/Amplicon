## Please atation:
## This word is about a function, and not a current codes.
## You can use this codes to get your results: `NCM(NCM(data = otu, filename1 = "./NCM.pdf", filename2 = "./NCM_pie.pdf"))`

NCM <- function(data, filename1, filename2) {
  library(minpack.lm)
  library(grid)
  library(ggplot2)
  data <- t(data)
  N <- mean(apply(data, 1, sum))
  p.m <- apply(data, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m / N
  
  data.bi <- 1 * (data > 0)
  freq <- apply(data.bi, 2, mean)
  freq <- freq[freq != 0]
  
  C <- merge(p, freq, by = 0)
  C <- C[order(C[, 2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  p <- C.0[, 2]
  freq <- C.0[, 3]
  names(p) <- C.0[, 1]
  names(freq) <- C.0[, 1]
  
  d <- 1 / N
  m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))
  m.ci <- confint(m.fit, 'm', level = 0.95)
  freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
  pred.ci <- binconf(freq.pred * nrow(data), nrow(data), alpha = 0.05, method = "wilson", return.df = TRUE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
  bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[, 2:3])
  bacnlsALL$group <- "no-significant"
  bacnlsALL$group[bacnlsALL$freq <= bacnlsALL$Lower] <- "Lower"
  bacnlsALL$group[bacnlsALL$freq >= bacnlsALL$Upper] <- "Higher"
  inter.col <- rep('black', nrow(bacnlsALL))
  inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- '#A52A2A'
  inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- '#29A6A6'
  
  pdf(file = filename1, width = 8, height = 6)
  grid.newpage()
  pushViewport(viewport(h = 0.6, w = 0.6))
  pushViewport(dataViewport(xData = range(log10(bacnlsALL$p)), yData = c(0, 1.02), extension = c(0.02, 0)))
  grid.rect()
  grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, gp = gpar(col = inter.col, cex = 0.7))
  grid.yaxis()
  grid.xaxis()
  grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp = gpar(col = 'blue', lwd = 2), default = 'native')
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp = gpar(col = 'blue', lwd = 2, lty = 2), default = 'native')
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp = gpar(col = 'blue', lwd = 2, lty = 2), default = 'native')
  grid.text(y = unit(0, 'npc') - unit(2.5, 'lines'), label = 'Mean Relative Abundance (log10)', gp = gpar(fontface = 2))
  grid.text(x = unit(0, 'npc') - unit(3, 'lines'), label = 'Frequency of Occurance', gp = gpar(fontface = 2), rot = 90)
  dev.off()
  
  group_counts <- table(bacnlsALL$group)
  group_data <- as.data.frame(group_counts)
  colnames(group_data) <- c("Group", "Count")
  cols <- c("Higher" = "#29A6A6", "Lower" = "#A52A5A", "no-significant" = "black")
  pie_plot <- ggplot(group_data, aes(x = "", y = Count, fill = Group)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(label = paste0(round(Count / sum(Count) * 100, 1), "%")), position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols) +
    theme_void() +
    labs(fill = "Group", title = "")
  ggsave(filename = filename2, plot = pie_plot, width = 3, height = 3)
  return(list(model = m.fit, Rsqr = Rsqr, bacnlsALL = bacnlsALL))
}