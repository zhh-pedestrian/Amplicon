# Set the encoding format to UTF-8.

### 1 Data Analysis-----------------------
data <- read.csv(file = "./0 otu_0.1%过滤.csv")
tax <- read.csv(file = "./0 tax_0.1%过滤.csv")

data_list <- list(data = data, tax = tax)
saveRDS(data_list, file = "./Code/Data/Figrue 2.Sourcedata.rds")

data_list <- read.csv(file = "./Data/Figure 2.Sourcedata.rds")
data <- list$data
tax <- list$tax


otu_GJ_GM <- data[, 1:108]
otu_GJ_GM_1 <- round(otu_GJ_GM)
conditions_GJ_GM <-  data.frame(conditions = factor(c(rep("GJ", 54), rep("GM", 54))))
rownames(conditions_GJ_GM) <- colnames(otu_GJ_GM)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = otu_GJ_GM_1,
  colData = conditions_GJ_GM,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast = c("conditions","GJ","GM")
res = results(dds, contrast)

baseMean_GJ <- rowMeans(counts(dds, normalized = TRUE)[, colData(dds)$conditions == "GJ"])
baseMean_GM <- rowMeans(counts(dds, normalized = TRUE)[, colData(dds)$conditions == "GM"])
res = cbind(baseMean_GJ, baseMean_GM, as.data.frame(res))
res$padj[is.na(res$padj)] <- 1
res = as.data.frame(res[order(res$pvalue), ])

tax <- read.csv(file = "0 tax_0.1%过滤.csv")
res_1 <- merge(res, tax[, c("ID", "phylum")], by.x = "row.names", by.y = "ID")

### 2 Plot------------------------------

rm(list = ls())
my_theme <- theme_bw()+
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    panel.background = element_rect(fill = NULL),
    axis.title.x = element_text(vjust = 0.5, size = 15, color = "black"),
    axis.title.y = element_text(vjust = 1, size = 15, color = "black"),
    axis.text.x = element_text(size = 12, angle = 0, color = "black", vjust = 0.5, hjust = 0.4),
    axis.text.y = element_text(size = 12, angle = 0, color = "black", vjust = 0.5, hjust = 0),
    axis.line = element_line(colour = "black", size = 0.8, lineend = "butt"),
    axis.line.x.top = element_line(colour = "black", size = 0.8, lineend = "butt"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(colour = "black", size = 0.8)
  )

data <- res_1
data$change <- ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >1, ifelse(data$log2FoldChange < 0, "depleted", "enriched"), "nosing")
data$pvalue <- -log10(data$pvalue)
data$label <- ""
data$label <- ifelse(data$pvalue > 10, data$Row.names, "")

phylum_counts <- data %>%
  group_by(phylum) %>%
  summarise(count = n()) %>%
  ungroup()
valid_phylum <- phylum_counts %>%
  filter(count >= 20)
filtered_data <- data %>%
  filter(phylum %in% valid_phylum$phylum)

cols <- c("#E3716E", "#ECA680", "#7AC7E2", 
          "#F7DF87", "#54BEAA", "#CFAFD4", 
          "#E4A6BD", "#8CCBEA", "#2983B1")
ggplot(filtered_data, aes(phylum, pvalue))+
  geom_jitter(aes(color = phylum, shape = change, size = abs(log2FoldChange)), width = 0.5)+
  geom_label_repel(aes(label = label, color = phylum), size = 4, 
                   box.padding = unit(0.5, "lines"), point.padding = unit(0.8, "lines"), 
                   segment.color = "grey50", show.legend = FALSE, max.overlaps = 10000) +
  scale_y_continuous(expand = c(0,0),limits = c(0,30))+#更改y轴范围
  scale_size_continuous(range = c(1,3.5))+#更改点的大小极值
  scale_shape_manual(values = c(6,17,16))+
  scale_size_continuous(limits = c(0,15), breaks = c(0,3,8,15))+
  geom_hline(yintercept = 10, size = 1, linetype = "dashed", color = "grey70")+
  scale_color_manual(values = cols)+
  labs(x = "", y = "-Log(Pvalue)")+
  my_theme
