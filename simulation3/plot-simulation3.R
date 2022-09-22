data <- read.csv("Result.csv")
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)
data$Score <- (data$Power-data$FWER+0.5-0.5*data$EN/120)
data <- rename(data, BF_E = phi.e, BF_P = phi.p)

p1 <- ggplot(data, aes(x=BF_P, y=BF_E))+
  geom_raster(aes(fill = EN), interpolate = TRUE)+
  scale_fill_gradientn(colours = colormap)
  # theme(panel.background=element_rect(fill="white", color = "black"))+

p2 <- ggplot(data, aes(x=BF_P, y=BF_E))+
  geom_raster(aes(fill = FWER), interpolate = TRUE)+
  scale_fill_gradientn(colours = colormap)

p3 <- ggplot(data, aes(x=BF_P, y=BF_E))+
  geom_raster(aes(fill = Power), interpolate = TRUE)+
  scale_fill_gradientn(colours = colormap)

p4 <- ggplot(data, aes(x=BF_P, y=BF_E))+
  geom_raster(aes(fill = Score), interpolate = TRUE)+
  scale_fill_gradientn(colours = colormap)
a <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)
ggsave("Figure 6.jpeg", plot = a, device = "jpeg", units = "in", width = 12, height = 8.5, dpi = 600)
