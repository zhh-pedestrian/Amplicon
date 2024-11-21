# Set the encoding format to UTF-8.

bacillus <- readRDS(file = "./Data/Figure 4.Sourcedata.rds")

bacillus_G <- bacillus %>%
  mutate(GJ = rowMeans(select(., 1:54)),
         GM = rowMeans(select(., 55:108)),
         GN = rowMeans(select(., 109:162))) %>% 
  select(GJ, GM, GN)

bacillus_R <- bacillus %>%
  mutate(rank0 = rowMeans(select(., contains("R0"))),
         rank1 = rowMeans(select(., contains("R1"))),
         rank3 = rowMeans(select(., contains("R3"))),
         rank5 = rowMeans(select(., contains("R5"))),
         rank7 = rowMeans(select(., contains("R7"))),
         rank9 = rowMeans(select(., contains("R9")))) %>% 
  select(rank0, rank1, rank3, rank5, rank7, rank9)

bacillus_G_scale <- scale(bacillus_G)
bacillus_G_scale_long <- melt(bacillus_G_scale, id.vars = "row.names")
colnames(bacillus_G_scale_long) <- c("OTU", "group", "value")

colors1 <- brewer.pal(9, "Set1")   # 9种颜色
colors2 <- brewer.pal(12, "Set3")  # 12种颜色
colors <- c(colors1, colors2[1:(18 - length(colors1))])

colors <- colorRampPalette(brewer.pal(9, "Blues"))(18)
colors <- colorRampPalette(c("#EE6A50", "#fff8dc"))(18)
ggplot(bacillus_G_scale_long, aes(x = group, y = OTU, color = OTU))+
  geom_point(aes(size = value), alpha = 0.9)+
  scale_color_manual(values = colors) +
  scale_size_continuous(breaks = c(0.1, 0.5, 2, 10), labels = c(0.1, 0.5, 2, 10), range = c(1, 10))+
  guides(color = FALSE)+
  theme_bw()+
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