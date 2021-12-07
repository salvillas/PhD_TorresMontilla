library(tidyverse)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(Hmisc)
library(ggsci)
library(gridExtra)
library(readr)
library(seqinr)
library(VennDiagram)
library(phylotools)
library(ggpubr)

###############################################################################################
theme_HPLC <- theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                    axis.title = element_text(face = "bold"),
                    #strip.text = element_text(face = "bold"),
                    panel.background = element_blank(),
                    legend.background = element_blank(),
                    legend.key = element_blank(),
                    panel.grid = element_blank(),
                    axis.text = element_text(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank()
)

folders <- grep("Berlin", list.dirs(), value = TRUE)

All_PSII <- data.frame(GFP_1= as.numeric(),
                       GFP_2= as.numeric(),
                       GFP_3= as.numeric(),
                       GFP_4= as.numeric(),
                       crtB_1= as.numeric(),
                       crtB_2= as.numeric(),
                       crtB_3= as.numeric(),
                       crtB_4= as.numeric(),
                       TimePoint= as.character(),
                       Sample= as.character())

for (folder in folders){
  for (file in grep("PSII.csv", list.files(folder), value = TRUE)){
    PSII <- read.csv(paste(folder, file, sep = "/"), header = TRUE, sep = ";")
    PSII <- PSII[nrow(PSII),grep("Y", names(PSII), value = TRUE)]
    names(PSII) <- c("GFP_1", "GFP_2", "GFP_3", "GFP_4", "crtB_1", "crtB_2", "crtB_3", "crtB_4")
    PSII$TimePoint <- str_sub(folder, -5,-1)
    PSII$Sample <- str_sub(file, -12, -10)
    All_PSII <- rbind(All_PSII, PSII)
  }
}

All_PSII <- All_PSII %>% gather(key = "Infiltration", value = "PSII", -TimePoint, -Sample) %>%
  separate(Infiltration, into= c("Construct", "TechRep"), sep= "_") %>%
  separate(Sample, into= c("Replicate", "Leaf"), sep = "-")


All_PSII$TimePoint <- as.factor(All_PSII$TimePoint)
All_PSII$Construct <- factor(All_PSII$Construct, levels = c("GFP", "crtB"))
All_PSII$Replicate <- as.factor(All_PSII$Replicate)
All_PSII$Leaf <- as.factor(All_PSII$Leaf)
All_PSII$TechRep <- as.factor(All_PSII$TechRep)


PAM_ggplot <- ggplot(All_PSII, aes(x=TimePoint, y= PSII, fill= Construct))+
  geom_boxplot(alpha= 0.6)+
  geom_jitter(aes(col= Construct), position = position_jitterdodge(), alpha= 0.5)+
  scale_fill_manual(values = c("GFP" = "#4ba88c",
                               "crtB" = "#d46b42"))+
  scale_color_manual(values = c("GFP" = "#66c2a5",
                                "crtB" = "#fc8d62"))+
  #scale_fill_manual(values = c("GFP"= NULL,
  #                             "crtB" = NULL))+
  #scale_color_manual(values = c("GFP"= "#9ce3dc",
  #                              "crtB" = "#fdd381"))+
  ylab(expression(phi*"PSII"))+
  xlab("Time Points")+
  #geom_smooth(aes(group = Construct, col= Construct), method = "lm", se= FALSE, linetype= "dashed", size= 1.5)+
  theme_HPLC+
  theme(axis.text.x = element_text( size = 10),
        axis.title.x= element_text(face = "bold", size = 12),
        axis.title.y= element_text(face = "bold", size = 12),
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey"),
        legend.title = element_text(face = "bold"))
PAM_ggplot

PAM_ggplot <- PAM_ggplot + stat_compare_means(aes(group= Construct), label = "p.signif")



