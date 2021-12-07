library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(Hmisc)
library(ggsci)
library(gridExtra)
library(readr)

#Para 

HPLC <- read_delim("Prueba para R.txt", delim = "\t", col_names = TRUE)

HPLC2 <- HPLC %>% gather(sample, values, -Carotenoids) %>% separate(sample, into = c("Construct", "Time_point", "Replicate"), convert = TRUE, sep = "_")
HPLC3 <- HPLC2 %>% spread(Carotenoids, values) 
View(HPLC3)
HPLC3


HPLC4 <- read_delim("Prueba2 R.txt", delim = "\t", col_names = TRUE)
HPLC5 <- HPLC4 %>% 
  gather(sample, values, -Carotenoids) %>% 
  separate(sample, into = c("Construct", "Time_point", "Replicate"), convert = TRUE, sep = "_") %>%
  spread(Carotenoids, values) 

HPLC6 <- HPLC5 %>% 
  select(Construct, Time_point, Phytoene) %>%
  group_by(Construct, Time_point) %>%
  summarise(Mean = mean(Phytoene), SD = sd(Phytoene)) %>%
  mutate(top = Mean + SD, bottom = Mean - SD)

Phytoene <- ggplot(HPLC5, aes(x = Time_point, y = Phytoene, fill = Construct))+
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", alpha = 0.5)+
  geom_point(aes(col = Construct),shape = 20, size = 3, position = position_dodge(width = 0.9))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", position = "dodge", fun.args = list(mult = 1), width = 0.4)+
  #geom_errorbar(data = (HPLC5 %>% 
                          #select(Construct, Time_point, Phytoene) %>%
                          #group_by(Construct, Time_point) %>%
                          #summarise(Mean = mean(Phytoene), SD = sd(Phytoene))%>%
                          #mutate(top = Mean + SD, bottom = Mean - SD)),
                #aes(x = Time_point, y = Mean, ymin = bottom, ymax = top), width = 0.5)+
  scale_y_continuous("µg/mg dry weight", expand = expand_scale(mult = .01))+
  scale_x_discrete(limits = c("0h", "6h", "24h", "48h", "72h","96h"))+
  facet_grid(. ~ Construct)+
  scale_fill_manual(values = c("(p)crtB" = "#afafaf", "GFP" = "#cdcdcd"))+
  labs(title = "Phytoene", x = "Time points")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())
  
Phytoene

Lutein <- ggplot(HPLC5, aes(x = Time_point, y = Lutein, fill = Construct))+
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", alpha = 0.5)+
  geom_point(aes(col = Construct),shape = 20, size = 3, position = position_dodge(width = 0.9))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", position = "dodge", fun.args = list(mult = 1), width = 0.4)+
  scale_y_continuous("µg/mg dry weight", expand = expand_scale(mult = .01))+
  scale_x_discrete(limits = c("0h", "6h", "24h", "48h", "72h","96h"))+
  facet_grid(. ~ Construct)+
  scale_fill_manual(values = c("(p)crtB" = "#afafaf", "GFP" = "#cdcdcd"))+
  labs(title = "Lutein", x = "Time points")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())

Lutein

Beta.Carotene <- ggplot(HPLC5, aes(x = Time_point, y = `Beta-carotene total`, fill = Construct))+
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", alpha = 0.5)+
  geom_point(aes(col = Construct),shape = 20, size = 3, position = position_dodge(width = 0.9))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", position = "dodge", fun.args = list(mult = 1), width = 0.4)+
  scale_y_continuous("µg/mg dry weight", expand = expand_scale(mult = .01))+
  scale_x_discrete(limits = c("0h", "6h", "24h", "48h", "72h","96h"))+
  facet_grid(. ~ Construct)+
  scale_fill_manual(values = c("(p)crtB" = "#afafaf", "GFP" = "#cdcdcd"))+
  labs(title = "Beta-Carotene", x = "Time points")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())

Beta.Carotene

Lycopene <- ggplot(HPLC5, aes(x = Time_point, y = Lycopene, fill = Construct))+
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", alpha = 0.5)+
  geom_point(aes(col = Construct),shape = 20, size = 3, position = position_dodge(width = 0.9))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", position = "dodge", fun.args = list(mult = 1), width = 0.4)+
  scale_y_continuous("µg/mg dry weight", expand = expand_scale(mult = .01))+
  scale_x_discrete(limits = c("0h", "6h", "24h", "48h", "72h","96h"))+
  facet_grid(. ~ Construct)+
  scale_fill_manual(values = c("(p)crtB" = "#afafaf", "GFP" = "#cdcdcd"))+
  labs(title = "Lycopene", x = "Time points")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())

Lycopene

Chlorophyll.a <- ggplot(HPLC5, aes(x = Time_point, y = `Chlorophyll a`, fill = Construct))+
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", alpha = 0.5)+
  geom_point(aes(col = Construct),shape = 20, size = 3, position = position_dodge(width = 0.9))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", position = "dodge", fun.args = list(mult = 1), width = 0.4)+
  scale_y_continuous("µg/mg dry weight", expand = expand_scale(mult = .01))+
  scale_x_discrete(limits = c("0h", "6h", "24h", "48h", "72h","96h"))+
  facet_grid(. ~ Construct)+
  scale_fill_manual(values = c("(p)crtB" = "#afafaf", "GFP" = "#cdcdcd"))+
  labs(title = "Chlorophyll a", x = "Time points")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())

Chlorophyll.a

Chlorophyll.b <- ggplot(HPLC5, aes(x = Time_point, y = `Chlorophyll b`, fill = Construct))+
  geom_bar(stat = "summary", fun.y = "mean", position = "dodge", alpha = 0.5)+
  geom_point(aes(col = Construct),shape = 20, size = 3, position = position_dodge(width = 0.9))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar", position = "dodge", fun.args = list(mult = 1), width = 0.4)+
  scale_y_continuous("µg/mg dry weight", expand = expand_scale(mult = .01))+
  scale_x_discrete(limits = c("0h", "6h", "24h", "48h", "72h","96h"))+
  facet_grid(. ~ Construct)+
  scale_fill_manual(values = c("(p)crtB" = "#afafaf", "GFP" = "#cdcdcd"))+
  labs(title = "Chlorophyll b", x = "Time points")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())

Chlorophyll.b
