library(stringr)
library(dplyr)
library(ggplot2)
library(magrittr)
library(ggrepel)

ProteinCox<-read.csv("~/Data.csv")
ProteinCox <- ProteinCox %>%
  mutate(exposure = str_extract(exposure, "^[^;]+"))

ProteinCox$logP<- -log(ProteinCox$Pvalue, 10)
ProteinCox$regulate<-'not Sig'
ProteinCox$regulate[ProteinCox$logP>= 1 & ProteinCox$exp_coef>= 1]<-'up'
ProteinCox$regulate[ProteinCox$logP>= 1 & ProteinCox$exp_coef< 1]<-'down'
head(ProteinCox)

p <- ggplot(ProteinCox, aes(x = exp_coef, y = logP)) +
  geom_point(data = ProteinCox[ProteinCox$regulate == "up", ],  
             aes(color = regulate), size = 5, alpha = 1) +  
  geom_point(data = ProteinCox[ProteinCox$regulate == "down", ], 
             aes(color = regulate), size =5 , alpha = 1) +  
  geom_point(data = ProteinCox[ProteinCox$regulate == "not Sig", ], 
             aes(color = regulate), size =5 , alpha = 1) +  
  scale_color_manual(values = c('#A0B0D1','grey','#D39FAA')) + 
  scale_x_continuous(breaks = seq(0, 6, 1), limits = c(-0.2, 6.5)) +  
  scale_y_continuous(breaks = seq(0, 90, 20), limits = c(0, 90)) +  
  geom_vline(xintercept = c(1), linetype = 'dashed',size=0.5,alpha = 1) +  
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed',size=0.5,alpha = 1)+ 
  geom_hline(yintercept = -5, linetype = 'solid', size=0.2,alpha = 1)+ 
  labs(x = "Effect size (HR)", y = expression(-log[10](italic(P) ~ value)), parse = TRUE,size = 7) +
  ggtitle("Incident NAFLD (model 1)") +
  theme(plot.title = element_text(size = 5, hjust = 0))+
  theme_light()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(p)

newdata1<-ProteinCox%>%filter(exp_coef>=1 & Pvalue<0.05)%>%top_n(6,logP)
newdata2<-ProteinCox%>%filter(exp_coef<1 & Pvalue<0.05)%>%top_n(5,logP)

p2 <- p + 
  geom_point(data = newdata1,aes(x = exp_coef, y = -log10(Pvalue)),
             color = '#B15661', size = 8, alpha = 0.2) +
  geom_text_repel(data = newdata1,aes(x = exp_coef, y = logP, label = exposure,family = "Candara"),
                  seed = 23456, color = '#873942', show.legend = FALSE, 
                  max.segment.length = unit(0.3, "lines"), 
                  min.segment.length = 0, 
                  segment.linetype = 1, 
                  force = 0.21,
                  force_pull = 3,
                  size = 5,
                  box.padding = unit(1, "lines"),
                  point.padding = unit(0.1, "lines"),
                  max.overlaps = Inf) 
p2 <- p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)


p3 <- p2 + 
  geom_point(data = newdata2,aes(x = exp_coef, y = -log10(Pvalue)),
             color = '#2D5685', size = 8, alpha = 0.2) +
  geom_text_repel(data = newdata2,aes(x = exp_coef, y = logP, label = exposure,family = "Candara"),
                  seed = 23456, color = '#2D5685', show.legend = FALSE, 
                  max.segment.length = unit(0.05, "lines"), 
                  min.segment.length = 0, 
                  segment.linetype = 1, 
                  force = 0.05,
                  force_pull = 20,
                  size = 5,
                  box.padding = unit(1, "lines"),
                  point.padding = unit(0.1, "lines"),
                  max.overlaps = Inf)
p3 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size = 15, hjust = 0, family = "Candara"),
        axis.title = element_text(size = 14,family = "Candara"), 
        axis.text = element_text(size = 12,color = 'black', family = "Candara"), 
        legend.title = element_text(size = 14,family = "Candara"), 
        legend.text = element_text(size = 14,family = "Candara")) 
p3 <- p3 + theme(
  panel.border = element_blank(),       
  panel.background = element_blank(),   
  axis.line = element_line(size = 0.25, color = "black"), 
  axis.ticks.length = unit(0.25, "cm")  
)
p3 <- p3 + theme(
  axis.ticks = element_line(color = "black")  
)
print(p3)

