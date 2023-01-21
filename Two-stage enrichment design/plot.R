data <- read.csv("Result.csv")
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)
data$Score <- (data$Power-data$FWER+0.5-0.5*data$EN/120)  ## Pre-specified score function
data <- rename(data, BF_E = phi.e, BF_P = phi.p)

p1 <- ggplot(data, aes(x=BF_P, y=BF_E))+
  geom_raster(aes(fill = EN), interpolate = TRUE)+
  scale_fill_gradientn(colours = colormap)+
  theme(axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))
  # theme(panel.background=element_rect(fill="white", color = "black"))+

p2 <- ggplot(data, aes(x=BF_P, y=BF_E))+
  geom_raster(aes(fill = FWER), interpolate = TRUE)+
  scale_fill_gradientn(colours = colormap)+
  theme(axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))

p3 <- ggplot(data, aes(x=BF_P, y=BF_E))+
  geom_raster(aes(fill = Power), interpolate = TRUE)+
  scale_fill_gradientn(colours = colormap)+
  theme(axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))

p4 <- ggplot(data, aes(x=BF_P, y=BF_E))+
  geom_raster(aes(fill = Score), interpolate = TRUE)+
  scale_fill_gradientn(colours = colormap)+
  theme(axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))
a <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)
ggsave("Figure 6.pdf", plot = a, device = "pdf", units = "mm", width = 170, height = 110, dpi = 300)
