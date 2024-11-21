# Set the encoding format to UTF-8.

### 1 Calculate Alpha Diversity
data_list <- readRDS(file = "./Code/Data/Figrue 1.Sourcedata.rds")
otu <- data_list$otu
group <- data_list$group
est <- estimateR(t(otu))
Richness <- est[1,]
Chao1  <- est[2, ]
ACE  <- est[4, ]
Shannon <- diversity(t(otu), index = "shannon")
Gini_simpson  <-vegan:: diversity(t(otu), index = 'simpson')
Simpson <- 1 - Gini_simpson
Goods_coverage <- 1 - rowSums(t(otu) == 1) / rowSums(t(otu))

alpha <- as.data.frame(cbind(Richness, Chao1, ACE, Shannon, Simpson, Goods_coverage))
Alpha_merge_group <- merge(alpha, group, by.x = "row.names", by.y = "sample")
colnames(Alpha_merge_group)[1] <- "Sample"

### 2 Plot--------------------------
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
    axis.ticks = element_line(colour = "black", size = 0.8),
    legend.position = "none"
  )

data <- Alpha_merge_group

color = c("#EA8379", "#7DAEE0", "#B395BD")
ggplot(data = data, aes(x = group1, y = Richness, color = group1))+
  geom_boxplot(position = "dodge", alpha = 0, outlier.size = 0.1, size = 0.5, width = 0.6)+
  stat_boxplot(geom = "errorbar", width = 0.25, size = 0.9, aes(color = group1))+
  geom_jitter(size = 2, width = 0.2, alpha = 0.5)+
  scale_color_manual(values = color)+
  scale_y_continuous(limits = c(0, 5500))+
  labs(x = "", y = "")+
  theme_bw() +
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
    axis.ticks = element_line(colour = "black", size = 0.8),
    legend.position = "none"
  )

