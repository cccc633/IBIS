library(ggplot2)
library(cowplot)
library(lemon)
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")
## MSE
MSE <- read.csv("MSE.csv")
MSE$Method <- factor(MSE$Method, c("IBIS", "Independent", "BHM"))
dataforplot<- MSE[MSE$n==10, ]
p1 <- ggplot(data=dataforplot,aes(x=Scenario, y=V1, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (1,1)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.45), breaks = seq(0,0.4,by=0.1))
p2 <- ggplot(data=dataforplot,aes(x=Scenario, y=V2, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (1,2)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.1,by=0.1))
p3 <- ggplot(data=dataforplot,aes(x=Scenario, y=V3, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (1,3)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.12), breaks = seq(0,0.1,by=0.1))
p4 <- ggplot(data=dataforplot,aes(x=Scenario, y=V4, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (1,4)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.12), breaks = seq(0,0.1,by=0.1))
p5 <- ggplot(data=dataforplot,aes(x=Scenario, y=V5, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (2,1)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.1,by=0.1))
p6 <- ggplot(data=dataforplot,aes(x=Scenario, y=V6, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (2,2)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.1,by=0.1))
p7 <- ggplot(data=dataforplot,aes(x=Scenario, y=V7, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (2,3)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.1,by=0.1))
p8 <- ggplot(data=dataforplot,aes(x=Scenario, y=V8, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (2,4)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.1,by=0.1))
p9 <- ggplot(data=dataforplot,aes(x=Scenario, y=V9, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (3,1)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.12), breaks = seq(0,0.1,by=0.1))
p10 <- ggplot(data=dataforplot,aes(x=Scenario, y=V10, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (3,2)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.12), breaks = seq(0,0.1,by=0.1))
p11 <- ggplot(data=dataforplot,aes(x=Scenario, y=V11, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (3,3)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.15), breaks = seq(0,0.1,by=0.1))
p12 <- ggplot(data=dataforplot,aes(x=Scenario, y=V12, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "MSE", title = "Subgroup (3,4)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0,0.45), breaks = seq(0,0.4,by=0.1))
a <- grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,ncol = 4, nrow = 3,position='top')
ggsave("Figure 3.pdf", plot = a, device = "pdf", units = "mm", width = 170, height = 110, dpi = 300)

## Bias
Bias <- read.csv("Bias.csv")
Bias$Method <- factor(Bias$Method, c("IBIS", "Independent", "BHM"))
dataforplot<- Bias[Bias$n==10, ]
p1 <- ggplot(data=dataforplot,aes(x=Scenario, y=V1, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (1,1)")+
  scale_x_continuous(breaks = seq(1,8))
p2 <- ggplot(data=dataforplot,aes(x=Scenario, y=V2, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (1,2)")+
  scale_x_continuous(breaks = seq(1,8))
p3 <- ggplot(data=dataforplot,aes(x=Scenario, y=V3, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (1,3)")+
  scale_x_continuous(breaks = seq(1,8))
p4 <- ggplot(data=dataforplot,aes(x=Scenario, y=V4, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (1,4)")+
  scale_x_continuous(breaks = seq(1,8))
p5 <- ggplot(data=dataforplot,aes(x=Scenario, y=V5, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (2,1)")+
  scale_x_continuous(breaks = seq(1,8))
p6 <- ggplot(data=dataforplot,aes(x=Scenario, y=V6, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (2,2)")+
  scale_x_continuous(breaks = seq(1,8))
p7 <- ggplot(data=dataforplot,aes(x=Scenario, y=V7, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (2,3)")+
  scale_x_continuous(breaks = seq(1,8))
p8 <- ggplot(data=dataforplot,aes(x=Scenario, y=V8, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (2,4)")+
  scale_x_continuous(breaks = seq(1,8))
p9 <- ggplot(data=dataforplot,aes(x=Scenario, y=V9, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (3,1)")+
  scale_x_continuous(breaks = seq(1,8))
p10 <- ggplot(data=dataforplot,aes(x=Scenario, y=V10, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (3,2)")+
  scale_x_continuous(breaks = seq(1,8))
p11 <- ggplot(data=dataforplot,aes(x=Scenario, y=V11, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (3,3)")+
  scale_x_continuous(breaks = seq(1,8))
p12 <- ggplot(data=dataforplot,aes(x=Scenario, y=V12, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "Bias", title = "Subgroup (3,4)")+
  scale_x_continuous(breaks = seq(1,8))
a <- grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6,p7, p8,p9, p10,p11, p12,ncol = 4, nrow = 3,position='top')
ggsave("Figure 4.pdf", plot = a, device = "pdf", units = "mm", width = 170, height = 110, dpi = 300)

## Width
Width <- read.csv("Width.csv")
Width$Method <- factor(Width$Method, c("IBIS", "Independent", "BHM"))
dataforplot<- Width[Width$n==10, ]
p1 <- ggplot(data=dataforplot,aes(x=Scenario, y=V1, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (1,1)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p2 <- ggplot(data=dataforplot,aes(x=Scenario, y=V2, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (1,2)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p3 <- ggplot(data=dataforplot,aes(x=Scenario, y=V3, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (1,3)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p4 <- ggplot(data=dataforplot,aes(x=Scenario, y=V4, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (1,4)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p5 <- ggplot(data=dataforplot,aes(x=Scenario, y=V5, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (2,1)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p6 <- ggplot(data=dataforplot,aes(x=Scenario, y=V6, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (2,2)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p7 <- ggplot(data=dataforplot,aes(x=Scenario, y=V7, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (2,3)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p8 <- ggplot(data=dataforplot,aes(x=Scenario, y=V8, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (2,4)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p9 <- ggplot(data=dataforplot,aes(x=Scenario, y=V9, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (3,1)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p10 <- ggplot(data=dataforplot,aes(x=Scenario, y=V10, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (3,2)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p11 <- ggplot(data=dataforplot,aes(x=Scenario, y=V11, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (3,3)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
p12 <- ggplot(data=dataforplot,aes(x=Scenario, y=V12, color=Method))+
  geom_line(size=0.4) + geom_point(size=1, aes(shape=Method))+
  theme(legend.position = 'none',
        panel.background=element_rect(fill="white", color = "black"),
        plot.title = element_text(size = 6, hjust = 0.5),
        axis.text=element_text(size = 5),
        axis.title=element_text(size = 6),
        legend.text=element_text(size = 6),
        legend.title=element_text(size = 6))+
  labs(x = "Scenario", y = "95% CI Width", title = "Subgroup (3,4)")+
  scale_x_continuous(breaks = seq(1,8))+scale_y_continuous(limits = c(0.5,1.5), breaks = seq(0.5,1.5,by=0.5))
a <- grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6,p7, p8,p9, p10,p11, p12,ncol = 4, nrow = 3,position='top')
ggsave("Figure 5.pdf", plot = a, device = "pdf", units = "mm", width = 170, height = 110, dpi = 300)