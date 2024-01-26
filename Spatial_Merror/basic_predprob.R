rm(list = ls())
library(ggplot2)
library(forcats)
library(RColorBrewer)

setwd("D:/Study new/PhD/R code/Spatial_Merror")
load("predicted_probs.Rdata")

my.palette <- brewer.pal(n = 20, name = "YlGn")

levels(all_probs$names)
all_probs$names <- factor(all_probs$names, levels = levels(all_probs$names)[c(2,3,5,8,7,1,4,6)])
levels(all_probs$names)


# ggplot(all_probs , aes(x=Longitude,y=Latitude,fill=Probability)) +
#   geom_raster() +
#   scale_fill_gradientn(colours = my.palette) +
#   facet_wrap(. ~ names, nrow = 2, labeller = label_parsed) +
#   theme_bw() +
#   coord_cartesian(xlim=c(-100,-65)) +
#   theme(
#     legend.position="bottom", 
#     #axis.title=element_text(size=25), 
#     strip.text = element_text(size=12),
#     #axis.text=element_text(size=25),aspect.ratio = 1,
#     #axis.ticks.length=unit(.25, "cm")) +
#     #legend.key.size = unit(2, 'cm'), 
#     #legend.key.height = unit(2, 'cm'), 
#     #legend.key.width = unit(2, 'cm'),
#     #legend.title = element_text(size=5), 
#     legend.text = element_text(size=10)
#     ) 
# 
# ggsave("predicted_prob.pdf", width = 28, height = 14)

ggplot(all_probs , aes(x=Longitude,y=Latitude,z=Probability)) +
  geom_contour_filled(breaks = seq(0,0.45,by=0.05)) +
  scale_fill_brewer(palette = "YlGn") +
  facet_wrap(. ~ names, nrow = 2, labeller = label_parsed) +
  theme_bw() +
  labs(fill = "Predicted \nprobability") +
  coord_cartesian(xlim=c(-100,-65)) +
  theme(
    legend.position="bottom", 
    #axis.title=element_text(size=25), 
    strip.text = element_text(size=12),
    #axis.text=element_text(size=25),aspect.ratio = 1,
    #axis.ticks.length=unit(.25, "cm")) +
    #legend.key.size = unit(2, 'cm'), 
    #legend.key.height = unit(2, 'cm'), 
    #legend.key.width = unit(2, 'cm'),
    #legend.title = element_text(size=5), 
    legend.text = element_text(size=10)
    ) 


#ggsave("predicted_probs_contour.pdf", width = 10, height = 7)

#NB I need the labels to be that big to get it looking right in the pdf/manuscript
#Yeah the trick I use here is to play the reverse uno card i.e., the smaller the image the more magnified the text.


#####-------------dFRK only--------#######

all_probs %>% filter(method == "Bcorrect") %>%
  ggplot(aes(x=Longitude,y=Latitude,z=Probability)) +
  geom_contour_filled(breaks = seq(0,0.45,by=0.05)) +
  scale_fill_brewer(palette = "YlGn") +
  #facet_wrap(. ~ names, nrow = 2, labeller = label_parsed) +
  theme_bw() +
  labs(fill = "Predicted \nprobability",title = "Double FRK Predicted Probabilities") +
  coord_cartesian(xlim=c(-100,-65)) +
  theme(
    legend.position="right", 
    #axis.title=element_text(size=25), 
    strip.text = element_text(size=12),
    #axis.text=element_text(size=25),aspect.ratio = 1,
    #axis.ticks.length=unit(.25, "cm")) +
    #legend.key.size = unit(2, 'cm'), 
    #legend.key.height = unit(2, 'cm'), 
    #legend.key.width = unit(2, 'cm'),
    #legend.title = element_text(size=5), 
    legend.text = element_text(size=10),
    plot.title = element_text(hjust = 0.5)
  ) 

ggsave("dFRK_pred.pdf",width = 7,height = 7)
