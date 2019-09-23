# figure 2 for anthro paper
library(tidyverse)
library(scales) # sweet break formatting functions
library(ggpubr)
library(ggthemes)

samples_ign_mean <- read.csv("samples_ign_mean_Adam.csv") %>%
  dplyr::select(-X) 

samples_ign_mean$ign <- ifelse(samples_ign_mean$ign == "Human",
                               "Primarily Anthropogenic (>75%)",
                               "Primarily Lightning (>75%)")

means_sds <- samples_ign_mean %>%
  group_by(ign) %>%
  summarise(mean_frp = median(log(Mean_FRP_MODIS_mean), na.rm=T),
            sd_frp = sd(log(Mean_FRP_MODIS_mean), na.rm = T),
            mean_area = median(log(Mean_area_Short_mean)),
            sd_area = sd(log(Mean_area_Short_mean)),
            mean_jd = median(log(Std_JD_Short_mean2), na.rm=T),
            sd_jd = sd(log(Std_JD_Short_mean2), na.rm = T),
            mean_n = median(log(Number_fires_Short_mean2), na.rm=T),
            sd_n = sd(log(Number_fires_Short_mean2), na.rm = T)) %>%
  ungroup %>%
  mutate(lowerf = mean_frp - sd_frp,
         upperf = mean_frp + sd_frp,
         lowera = mean_area - sd_area,
         uppera = mean_area + sd_area,
         lowerj = mean_jd - sd_jd,
         upperj = mean_jd + sd_jd,
         lowern = mean_n - sd_n,
         uppern = mean_n + sd_n)

p1 <- ggplot(samples_ign_mean, aes(x=(Mean_area_Short_mean), 
                             y=(Mean_FRP_MODIS_mean),
                             color = ign)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_errorbar(data = means_sds, width=0, #alpha = .2,
                aes(x = exp(mean_area),
                    ymin = exp(mean_frp - sd_frp),
                    ymax = exp(mean_frp + sd_frp)),
                inherit.aes = FALSE) +
  geom_errorbarh(data = means_sds, height=0, #alpha = .2,
                 aes(y = exp(mean_frp),
                    xmin = exp(lowera),
                    xmax = exp(mean_area + sd_area)),
                inherit.aes = FALSE)+
  geom_point(data=means_sds, alpha = .2, size=1.5, aes(x=exp(mean_area), y=exp(mean_frp)), inherit.aes = FALSE) +
   ggtitle("a.")      +
   scale_y_continuous(trans = "log10") +
   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"),name = "Ignition Source")+
  ylab("Average Intensity (MW)") +
  xlab("Average Fire Size (ha)") +
  #annotation_logticks()  +
  theme_pubr();p1


p2 <- ggplot(samples_ign_mean, aes(x=Number_fires_Short_mean2, 
                     y=Std_JD_Short_mean2,
                     color = ign)) +
  geom_point(alpha = 0.4, size = 1) +
 
  geom_errorbar(data = means_sds, width=0, #alpha = .2,
                aes(x = exp(mean_n),
                    ymin = exp(mean_jd - sd_jd),
                    ymax = exp(mean_jd + sd_jd)),
                inherit.aes = FALSE) +
  geom_errorbarh(data = means_sds, height=0, #alpha = .2,
                 aes(y = exp(mean_jd),
                     xmin = exp(lowern),
                     xmax = exp(mean_n + sd_n)),
                 inherit.aes = FALSE)+
  ggtitle("b.")              + 
  geom_point(data=means_sds, alpha = .2, size=1.5, aes(x=exp(mean_n), y=exp(mean_jd)), inherit.aes = FALSE) +
  scale_y_continuous(trans = "log10") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +   scale_color_manual(values = c("#E69F00", "#56B4E9"),name = "Ignition Source")+
  xlab("Fire frequency (n fires)") +
  ylab("Season Length (days)") +
  #annotation_logticks()  +
  theme_pubr();p2

ggarrange(p1,p2, ncol=2, nrow=1, legend=c("right"), common.legend=TRUE)

