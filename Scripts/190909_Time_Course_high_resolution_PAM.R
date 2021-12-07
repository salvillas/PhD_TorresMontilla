library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(Hmisc)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(tidyverse)

mypal = pal_npg("nrc", alpha = 0.7)(10)
mypal = append(mypal,"#B3A86A")

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}



##### Since I cannot find the HPLC plot from Cornell data, I redo it on 21-11-2021

HPLC_Cornell <- read.table("190909 Time-Course high resolution PAM R.txt",
                           header = TRUE, sep = "\t", dec = ",")
HPLC_Cornell_ggplot <- HPLC_Cornell %>% gather(key = "Sample", value = "values", -Metabolite) %>%
  separate(Sample, into = c("Construct", "TimePoint", "Replicate"), sep = "_")

HPLC_Cornell_ggplot$Construct <- factor(HPLC_Cornell_ggplot$Construct, 
                                        levels = c("GFP", "PcrtB"),
                                        labels = c("GFP", "(p)crtB"))
HPLC_Cornell_ggplot$TimePoint <- factor(HPLC_Cornell_ggplot$TimePoint, labels = paste0(unique(HPLC_Cornell_ggplot$TimePoint), "pi"))
HPLC_Cornell_ggplot$Replicate <- as.factor(HPLC_Cornell_ggplot$Replicate)
HPLC_Cornell_ggplot <- HPLC_Cornell_ggplot %>% filter(Metabolite %in% c("Phytoene", "cisLycopene", "Betacarotene", "Lutein", "TotalViolaxanthin", "Neoxanthin", "TotalCarotenoids", "TotalTocopherols"))
HPLC_Cornell_ggplot$Metabolite <- factor(HPLC_Cornell_ggplot$Metabolite,
                                         levels = c("Phytoene", "cisLycopene", "Betacarotene", "Lutein", "TotalViolaxanthin", "Neoxanthin", "TotalCarotenoids", "TotalTocopherols"),
                                         labels = c("Phytoene", "cis-Lycopene", "Betacarotene", "Lutein", "Violaxanthin", "Neoxanthin", "Total Carotenoids", "Total Tocopherols"))

All <- ggplot(HPLC_Cornell_ggplot, aes(x = TimePoint, y = values, fill = Construct))+
  geom_bar(stat = "summary", fun.y = "mean", position = position_dodge())+
  geom_point(shape = 20, size = 3, aes(col = Construct), position = position_dodge(width = 0.9))+
  stat_summary(fun.data = mean_sdl, geom = "errorbar",fun.args = list(mult = 1), width = 0.4, position = position_dodge(width=0.9))+
  ylab(expression(mu*"g/mg dry weight"))+
  scale_y_continuous(expand = expand_scale(mult = c(0, .2)))+
  #scale_x_discrete(labels = c("GFP", "(c)GFP", "1S-crtB", "crtB", "(p)crtB"))+
  scale_fill_manual(values = c("GFP"= "#66c2a5",
                               "(p)crtB" = "#fc8d62"))+
  scale_color_manual(values = c("GFP"= "#1b9e77",
                                "(p)crtB" = "#d95f02"))+
  labs(fill = "Constructs", col = "Constructs")+
  geom_line(stat = "smooth",aes(group = Construct, col= Construct), method = "lm", se= FALSE, linetype= "dotted", size= 1.5, alpha= 0.8)+
  facet_wrap(~Metabolite, scales = "free_y")+
  theme_HPLC+
  theme(axis.text.x = element_text(face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 8),
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey"))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
All





#################################################################################


HPLC1 <- read.delim("190909 Time Course high resolution PAM.txt",
                    header = TRUE, dec = ".", sep ="\t")
HPLC1$Replicate <- as.factor(HPLC1$Replicate)
str(HPLC1)

PAM <- read.delim("Time_course_PAM_ggplot2.txt", 
                  header = TRUE, dec = ".", sep = "\t")
PAM$mean <- ave(PAM$PSII, PAM$Replicate, PAM$Time_Point, PAM$Construct)

PAM <- PAM %>%
  select(Construct, Time_Point, Replicate, mean)
PAM <- unique(PAM)

HPLC2 <- merge(HPLC1, PAM)
colnames(HPLC2)[15] <- "PAM"

HPLC3 <- HPLC2 %>%
  filter(Construct=="crtB")
View(HPLC3)
RatioGFP_crtB <- c(1.039490073,	1.02740982,	1.108775335,	1.113230427,	1.037218443,	1.160039279,	1.259112194,	1.295319576,	1.237263324,	1.343875795,	1.27281035,	1.545544281,	1.609301209,	1.532089881,	1.800601081,	1.646510866,	1.60535929,	1.43645169,	1.552836012,	1.956059745,	1.911111675,	0,	1.89573529,	1.708533399,	1.974668404,	1.836579918,	1.846385554,	2.003859587,	2.086141713,	2.309712041,	2.122031268,	2.166774508,	2.021797233)
HPLC3$RatioTotCar <- RatioGFP_crtB

PAM_Ratio_Ave <- PAM_Ratio %>%
  group_by(Time_Point, Replicate) %>%
  summarise(AvePSII= mean(Ratio))
PAM_Ratio_Ave$TP <- as.numeric(str_sub(PAM_Ratio_Ave$Time_Point, 1, -2))
PAM_Ratio_Ave$Time_Point <- paste0(PAM_Ratio_Ave$TP-2, "hpi")
PAM_Ratio_Ave %>% select(-TP)
PAM_Ratio_Ave$ID <- paste(PAM_Ratio_Ave$Time_Point, PAM_Ratio_Ave$Replicate, sep = "_")

HPLC_Ratio <- HPLC1 %>% gather(key = "Metabolite", value = "values", -Construct, -Time_Point, -Replicate) %>%
  spread(key = Construct, value = values)  %>%
  mutate(Ratio = crtB/GFP) %>%
  select(Time_Point, Replicate, Metabolite, Ratio) %>%
  spread(key = Metabolite, value = Ratio)
HPLC_Ratio$Time_Point <- paste0(HPLC_Ratio$Time_Point, "pi")
HPLC_Ratio$ID <- paste(HPLC_Ratio$Time_Point, HPLC_Ratio$Replicate, sep = "_")

HPLC_PAM_Ratio <- merge(HPLC_Ratio, PAM_Ratio_Ave[,c("ID", "AvePSII")], by= "ID", no.dups = TRUE)
HPLC_PAM_Ratio <- HPLC_PAM_Ratio %>% select(-ID)

phytoene <- ggplot(HPLC1, aes(x = Time_Point, y = Phytoene, shape = Construct, col = Replicate)) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_npg()+
  theme_minimal()+
  facet_grid(.~Construct)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "Grey"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 8))+
  labs(x = "Time Points", y="?g/mg dry weight", title = "Phytoene")
phytoene

lycopene <- ggplot(HPLC1, aes(x = Time_Point, y = Lycopene, shape = Construction, col = Replicate)) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_npg()+
  theme_minimal()+
  facet_grid(.~Construction)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "Grey"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 8))+
  labs(x = "Time Points", y="?g/mg dry weight", title = "Lycopene")
lycopene


lutein <- ggplot(HPLC1, aes(x = Time_Point, y = Lutein, shape = Construction, col = Replicate)) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_npg()+
  theme_minimal()+
  facet_grid(.~Construction)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "Grey"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 8))+
  labs(x = "Time Points", y="?g/mg dry weight", title = "Lutein")
lutein

carotene <- ggplot(HPLC1, aes(x = Time_Point, y = Beta.Carotene, shape = Construction, col = Replicate)) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_npg()+
  theme_minimal()+
  facet_grid(.~Construction)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "Grey"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 8))+
  labs(x = "Time Points", y="?g/mg dry weight", title = "Beta-Carotene")
carotene

totalcar <- ggplot(HPLC1, aes(x = Time_Point, y = Total_Carotenoids, shape = Construction, col = Replicate)) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_npg()+
  theme_minimal()+
  facet_grid(.~Construction)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "Grey"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 8))+
  labs(x = "Time Points", y="?g/mg dry weight", title = "Total Carotenoids")
totalcar

viola <-ggplot(HPLC1, aes(x = Time_Point, y = Violaxanthin, shape = Construction, col = Replicate)) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_npg()+
  theme_minimal()+
  facet_grid(.~Construction)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "Grey"),
        axis.title = element_text(face = "bold"))+
  labs(x = "Time Point", y="?g/mg dry weight", title = "Violaxanthin")
viola

neo<- ggplot(HPLC1, aes(x = Time_Point, y = Neoxanthin, shape = Construction, col = Replicate)) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_npg()+
  theme_minimal()+
  facet_grid(.~Construction)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "Grey"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 8))+
  labs(x = "Time Point", y="?g/mg dry weight", title = "Neoxanthin")
neo

totaltoc <- ggplot(HPLC1, aes(x = Time_Point, y = Total_Tocopherols, shape = Construction, col = Replicate)) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_npg()+
  theme_minimal()+
  facet_grid(.~Construction)+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "Grey"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 8))+
  labs(x = "Time Points", y="?g/mg dry weight", title = "Total Tocopherols")
totaltoc

phy_PAM <- ggplot(HPLC3, aes(x = PAM, y = Phytoene, col = Time_Point, shape = Replicate))+
  geom_point(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_manual(values = mypal, label= paste0(unique(HPLC3$Time_Point), "pi"))+
  scale_x_reverse()+
  scale_y_continuous(expand = expansion(mult = c(0, .2)))+
  labs(x = "Ratio PSII (p)crtB/GFP", y=expression("Phytoene ("*mu*"g/mg dry weight)"), title = "Phytoene vs PAM", col = "Time Point")+
  theme_HPLC+
  theme(axis.title.x = element_text(face = "bold"),
        title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 8),
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey"),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"))
phy_PAM
str(HPLC3)

totalcar_PAM <- ggplot(HPLC3, aes(x = PAM, y = Total_Carotenoids, col = Time_Point, shape = Replicate))+
  geom_point(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_manual(values = mypal)+
  scale_x_reverse()+
  scale_y_continuous(expand = expansion(mult = c(0, .2)))+
  labs(x = "Ratio PSII (p)crtB/GFP", y=expression("Total Carotenoids ("*mu*"g/mg dry weight)"), title = "Total Carotenoids vs PAM", col = "Time Point")+
  theme_HPLC+
  theme(axis.title.x = element_text(face = "bold"),
        title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 8),
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey"),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"))
totalcar_PAM

ratiototcar_PAM <- ggplot(HPLC_PAM_Ratio, aes(x = AvePSII, y = Total_Carotenoids, col = Time_Point, shape = Replicate))+
  geom_point(size = 3, width = 0.2, alpha = 0.8) +
  scale_colour_manual(values = mypal)+
  scale_x_reverse()+
  scale_y_continuous(expand = expansion(mult = c(0, .2)))+
  labs(x = expression("Ratio "*phi*"PSII (p)crtB/GFP"), y="Ratio Total Carotenoids (p)crtB/GFP",  col = "Time Point")+
  theme_HPLC+
  theme(axis.title.x = element_text(face = "bold"),
        title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 8),
        strip.text = element_text(size = 8),
        legend.position = "none",
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey"),
        panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"))
ratiototcar_PAM


##Datos relevantes para el Clustering del RNAseq en Cornell

HPLC_PAM_Cornell <- HPLC2 %>%
  select(Construct, Time_Point, Replicate, Phytoene, Total_Carotenoids, PAM)

write.table(HPLC_PAM_Cornell, "G:/Mi unidad/CRAG/RNAseq/3er RNAseq/RNAseq_Cornell_HPLC_PAM.txt",
            sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


theme_HPLC <- theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                    axis.title = element_text(face = "bold"),
                    strip.text = element_text(face = "bold"),
                    panel.background = element_blank(),
                    legend.background = element_blank(),
                    legend.key = element_blank(),
                    panel.grid = element_blank(),
                    axis.text = element_text(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank()
)

#####################################################################################
######################## PAM Boxplot for Presentation ###############################
setwd("G:/Mi unidad/CRAG/HPLC/190906 Time Course high resolution PAM")

PAM <- read.delim("/Users/salva/Google Drive/My Drive/CRAG/PAM/190731 Time Course/Time_course_PAM_ggplot2.txt", 
                  header = TRUE, dec = ".", sep = "\t")
PAM$Construct <- factor(PAM$Construct, levels= c("GFP", "(p)crtB"))
PAM$TecRep <- c(1,2,3)

PAM_ggplot <- ggplot(PAM, aes(x=Time_Point, y= PSII, fill= Construct))+
  geom_boxplot()+
  geom_jitter(aes(col= Construct), position = position_jitterdodge(), alpha= 0.5)+
  #scale_fill_manual(values = c("GFP" = "#1a9850",
  #                             "(p)crtB" = "#f46d43"))+
  #scale_color_manual(values = c("GFP" = "#66bd63",
  #                              "(p)crtB" = "#fdae61"))+
  scale_fill_manual(values = c("GFP"= "#66c2a5",
                               "(p)crtB" = "#fc8d62"))+
  scale_color_manual(values = c("GFP"= "#1b9e77",
                                "(p)crtB" = "#d95f02"))+
  ylab(expression(phi*"PSII"))+
  xlab("Time Points")+
  geom_line(stat = "smooth",aes(group = Construct, col= Construct), method = "lm", se= FALSE, linetype= "dashed", size= 1.5, alpha= 0.8)+
  theme_HPLC+
  theme(axis.text.x = element_text(face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 8),
        strip.text = element_text(size = 8),
        legend.position = "none",
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey"))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
PAM_ggplot

PAM_Ratio <- PAM %>% spread(key = Construct, value = PSII) %>%
  mutate(Ratio= `(p)crtB`/GFP)
PAM_Ratio$Replicate <- as.factor(PAM_Ratio$Replicate)

PAM_ggplot_2 <- ggplot(PAM_Ratio, aes(x=Time_Point, y= Ratio))+
  geom_jitter(aes(col= Replicate), position = position_jitterdodge(),size= 3)+
  #scale_fill_manual(values = c("GFP" = "#1a9850",
  #                             "(p)crtB" = "#f46d43"))+
  #scale_color_manual(values = c("GFP" = "#66bd63",
  #                              "(p)crtB" = "#fdae61"))+
  scale_colour_manual(values = mypal[1:3])+
  scale_x_discrete(label= paste0(c("22", "25", "28", "31", "34", "37", "40", "43", "46", "56", "68"), "hpi"))+
  ylab(expression("Ratio "*phi*"PSII"))+
  xlab("Time Points")+
  #geom_line(stat = "smooth",aes(group = Construct, col= Construct), method = "lm", se= FALSE, linetype= "dashed", size= 1.5, alpha= 0.8)+
  theme_HPLC+
  theme(axis.text.x = element_text(face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 8),
        strip.text = element_text(size = 8),
        legend.position = "none",
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey"))
PAM_ggplot_2

PAM_ggplot <- PAM_ggplot + stat_compare_means(aes(group= Construct), label = "p.signif")


###############################################################################################
theme_HPLC <- theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                    axis.title = element_text(face = "bold"),
                    strip.text = element_text(face = "bold"),
                    panel.background = element_blank(),
                    legend.background = element_blank(),
                    legend.key = element_blank(),
                    panel.grid = element_blank(),
                    axis.text = element_text(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    axis.ticks.x = element_blank()
)

#######################################################################################
HPLC_Cornell_crtB <- HPLC_PAM_Cornell[HPLC_PAM_Cornell$Construct=="crtB",] %>%
  select(Time_Point, Replicate, Phytoene, Total_Carotenoids)
HPLC_Cornell_crtB <- HPLC_Cornell_crtB[!(HPLC_Cornell_crtB$Time_Point=="43h" & HPLC_Cornell_crtB$Replicate=="1"),]

colnames(HPLC_Cornell_crtB) <- c("Time_Point", "Replicate", "Phytoene_crtB", "Carotenoids_crtB")


HPLC_Cornell_GFP <- HPLC_PAM_Cornell[HPLC_PAM_Cornell$Construct=="GFP",] %>%
  select(Time_Point, Replicate, Phytoene, Total_Carotenoids)
colnames(HPLC_Cornell_GFP) <- c("Time_Point", "Replicate", "Phytoene_GFP", "Carotenoids_GFP")

HPLC_Cornell_Ratio <- merge(HPLC_Cornell_crtB, HPLC_Cornell_GFP, by= c("Time_Point", "Replicate"))

HPLC_Cornell_Ratio$Phytoene_Ratio <- HPLC_Cornell_Ratio$Phytoene_crtB / HPLC_Cornell_Ratio$Phytoene_GFP
HPLC_Cornell_Ratio$Carotenoids_Ratio <- HPLC_Cornell_Ratio$Carotenoids_crtB / HPLC_Cornell_Ratio$Carotenoids_GFP

HPLC_Cornell_Ratio <- HPLC_Cornell_Ratio %>%
  select(Time_Point, Replicate, Phytoene_Ratio, Carotenoids_Ratio)

HPLC_Cornell_Ratio <- HPLC_Cornell_Ratio[!(HPLC_Cornell_Ratio$Time_Point=="31h" | HPLC_Cornell_Ratio$Time_Point=="43h" | HPLC_Cornell_Ratio$Time_Point=="68h"),]
HPLC_Cornell_Ratio <- HPLC_Cornell_Ratio[!(HPLC_Cornell_Ratio$Time_Point=="28h" & HPLC_Cornell_Ratio$Replicate=="2"),]
HPLC_Cornell_Ratio <- HPLC_Cornell_Ratio[!(HPLC_Cornell_Ratio$Time_Point=="34h" & HPLC_Cornell_Ratio$Replicate=="1"),]
HPLC_Cornell_Ratio <- HPLC_Cornell_Ratio[!(HPLC_Cornell_Ratio$Time_Point=="37h" & HPLC_Cornell_Ratio$Replicate=="2"),]

#####

HPLC_Briardo_96h <- read.delim("G:/Mi unidad/CRAG/HPLC/RNAseq Briardo/171214 BRIARDO N. benthamiana for RNAseq 96h.txt",sep = "\t", dec = ",", header = TRUE)
HPLC_Briardo_96h <- HPLC_Briardo_96h %>% 
  gather(sample, values, -Metabolite) %>% 
  separate(sample, into = c("Construct","TimePoint", "Replicate"), convert = TRUE, sep = "_") %>%
  spread(Metabolite, values)
HPLC_Briardo_96h <- HPLC_Briardo_96h %>% select(Construct, TimePoint, Replicate, Phytoene, Total_Carotenoids)

HPLC_Briardo_96h_crtB <- HPLC_Briardo_96h[HPLC_Briardo_96h$Construct=="crtB",] %>%
  select(TimePoint, Replicate, Phytoene, Total_Carotenoids)

colnames(HPLC_Briardo_96h_crtB) <- c("Time_Point", "Replicate", "Phytoene_crtB", "Carotenoids_crtB")

HPLC_Briardo_96h_GFP <- HPLC_Briardo_96h[HPLC_Briardo_96h$Construct=="GFP",] %>%
  select(TimePoint, Replicate, Phytoene, Total_Carotenoids)

colnames(HPLC_Briardo_96h_GFP) <- c("Time_Point", "Replicate", "Phytoene_GFP", "Carotenoids_GFP")

HPLC_Briardo_96h_Ratio <- merge(HPLC_Briardo_96h_crtB, HPLC_Briardo_96h_GFP, by= c("Time_Point", "Replicate"), all = TRUE)

HPLC_Briardo_96h_Ratio$Phytoene_Ratio <- HPLC_Briardo_96h_Ratio$Phytoene_crtB / HPLC_Briardo_96h_Ratio$Phytoene_GFP
HPLC_Briardo_96h_Ratio$Carotenoids_Ratio <- HPLC_Briardo_96h_Ratio$Carotenoids_crtB / HPLC_Briardo_96h_Ratio$Carotenoids_GFP

HPLC_Briardo_96h_Ratio <- HPLC_Briardo_96h_Ratio %>%
  select(Time_Point, Replicate, Phytoene_Ratio, Carotenoids_Ratio)

### For CRAG Seminar
HPLC_Cornell_Briardo_Ratio <- rbind(HPLC_Cornell_Ratio, HPLC_Briardo_96h_Ratio)

HPLC_Cornell_Briardo_Ratio$Time_Point <- factor(HPLC_Cornell_Briardo_Ratio$Time_Point, 
                                                levels = rev(c("22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h")),
                                                labels = rev(c("1st", "2nd", "2nd", "3rd", "3rd", "4th", "5th", "6th", "7th")))
str(HPLC_Cornell_Briardo_Ratio)

Phytoene <- ggplot(HPLC_Cornell_Briardo_Ratio, aes(x=Time_Point, y=Phytoene_Ratio, fill = Time_Point))+
  geom_bar(stat = "summary", fun.y = "mean", position = position_dodge(),  width = 0.45)+
  stat_summary(fun.data = mean_sdl, geom = "errorbar",fun.args = list(mult = 1), width = 0.2, position = position_dodge(width=0.9), alpha=0.5)+
  coord_flip()+
  ylab("Phytoene (p)crtB / GFP")+
  scale_fill_manual(values= c("1st" = '#59cac0',
                      "2nd" = '#87c8a8',
                      "3rd" = '#a6c58f', 
                      "4th" = '#bec376',
                      "5th" = '#d3c05b',
                      "6th" = '#e4bd3c',
                      "7th" = '#f4ba00'))+
  theme_HPLC+
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(color="#272727"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "#626262"),
        legend.position = "none")
Phytoene
## Saved with 4x3 inches as PDF


Carotenoids <- ggplot(HPLC_Cornell_Briardo_Ratio, aes(x=Time_Point, y=Carotenoids_Ratio, fill=Time_Point))+
  geom_bar(stat = "summary", fun.y = "mean", position = position_dodge(),  width = 0.45)+
  stat_summary(fun.data = mean_sdl, geom = "errorbar",fun.args = list(mult = 1), width = 0.2, position = position_dodge(width=0.9), alpha=0.5)+
  coord_flip()+
  ylab("Total Carotenoids (p)crtB / GFP")+
  scale_fill_manual(values= c("1st" = '#59cac0',
                              "2nd" = '#87c8a8',
                              "3rd" = '#a6c58f', 
                              "4th" = '#bec376',
                              "5th" = '#d3c05b',
                              "6th" = '#e4bd3c',
                              "7th" = '#f4ba00'))+
  theme_HPLC+
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(color="#272727"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "#626262"),
        legend.position = "none")
Carotenoids
## Saved with 4x3 inches as PDF

######################################################################
### For 20/01/2021 Labmeeting
###################
## All Time-Points
HPLC_Cornell_Briardo_Ratio <- rbind(HPLC_Cornell_Ratio, HPLC_Briardo_96h_Ratio)

HPLC_Cornell_Briardo_Ratio$Time_Point <- factor(HPLC_Cornell_Briardo_Ratio$Time_Point, 
                                                levels = rev(c("22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h")))
str(HPLC_Cornell_Briardo_Ratio)

Phytoene <- ggplot(HPLC_Cornell_Briardo_Ratio, aes(x=Time_Point, y=Phytoene_Ratio, fill = Time_Point))+
  geom_bar(stat = "summary", fun.y = "mean", position = position_dodge(),  width = 0.45)+
  stat_summary(fun.data = mean_sdl, geom = "errorbar",fun.args = list(mult = 1), width = 0.2, position = position_dodge(width=0.9), alpha=0.5)+
  coord_flip()+
  ylab("Phytoene (p)crtB / GFP")+
  scale_fill_manual(values = c("22h" = '#59cac0',
                               "25h" = '#7ec8ae',
                               "28h" = '#98c69b',
                               "34h" = '#adc589',
                               "37h" = '#bec376',
                               "40h" = '#cec162',
                               "46h" = '#dcbf4c',
                               "56h" = '#e8bc33',
                               "96h" = '#f4ba00'
                               ))+
  theme_HPLC+
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(color="#272727"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "#626262"),
        legend.position = "none")
Phytoene
## Saved with 4x3 inches as PDF


Carotenoids <- ggplot(HPLC_Cornell_Briardo_Ratio, aes(x=Time_Point, y=Carotenoids_Ratio, fill=Time_Point))+
  geom_bar(stat = "summary", fun.y = "mean", position = position_dodge(),  width = 0.45)+
  stat_summary(fun.data = mean_sdl, geom = "errorbar",fun.args = list(mult = 1), width = 0.2, position = position_dodge(width=0.9), alpha=0.5)+
  coord_flip()+
  ylab("Total Carotenoids (p)crtB / GFP")+
  scale_fill_manual(values= c("22h" = '#59cac0',
                              "25h" = '#7ec8ae',
                              "28h" = '#98c69b',
                              "34h" = '#adc589',
                              "37h" = '#bec376',
                              "40h" = '#cec162',
                              "46h" = '#dcbf4c',
                              "56h" = '#e8bc33',
                              "96h" = '#f4ba00'
  ))+
  theme_HPLC+
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(color="#272727"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "#626262"),
        legend.position = "none")
Carotenoids
## Saved with 4x3 inches as PDF

###################################
## Distance Matrix Groups
HPLC_Cornell_Briardo_Ratio <- rbind(HPLC_Cornell_Ratio, HPLC_Briardo_96h_Ratio)

HPLC_Cornell_Briardo_Ratio$Time_Point <- factor(HPLC_Cornell_Briardo_Ratio$Time_Point, 
                                                levels = rev(c("22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h")),
                                                labels = rev(c("22h", "YY", "ZZ", "40h", "3rd", "4th", "5th", "6th", "7th")))
str(HPLC_Cornell_Briardo_Ratio)

Phytoene <- ggplot(HPLC_Cornell_Briardo_Ratio, aes(x=Time_Point, y=Phytoene_Ratio, fill = Time_Point))+
  geom_bar(stat = "summary", fun.y = "mean", position = position_dodge(),  width = 0.45)+
  stat_summary(fun.data = mean_sdl, geom = "errorbar",fun.args = list(mult = 1), width = 0.2, position = position_dodge(width=0.9), alpha=0.5)+
  coord_flip()+
  ylab("Phytoene (p)crtB / GFP")+
  #scale_fill_manual(values= c("1st" = '#59cac0',
  #                    "2nd" = '#87c8a8',
  #                    "3rd" = '#a6c58f', 
  #                    "4th" = '#bec376',
  #                    "5th" = '#d3c05b',
  #                    "6th" = '#e4bd3c',
  #                    "7th" = '#f4ba00'))+
  theme_HPLC+
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(color="#272727"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "#626262"),
        legend.position = "none")
Phytoene
## Saved with 4x3 inches as PDF


Carotenoids <- ggplot(HPLC_Cornell_Briardo_Ratio, aes(x=Time_Point, y=Carotenoids_Ratio, fill=Time_Point))+
  geom_bar(stat = "summary", fun.y = "mean", position = position_dodge(),  width = 0.45)+
  stat_summary(fun.data = mean_sdl, geom = "errorbar",fun.args = list(mult = 1), width = 0.2, position = position_dodge(width=0.9), alpha=0.5)+
  coord_flip()+
  ylab("Total Carotenoids (p)crtB / GFP")+
  scale_fill_manual(values= c("1st" = '#59cac0',
                              "2nd" = '#87c8a8',
                              "3rd" = '#a6c58f', 
                              "4th" = '#bec376',
                              "5th" = '#d3c05b',
                              "6th" = '#e4bd3c',
                              "7th" = '#f4ba00'))+
  theme_HPLC+
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(color="#272727"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "#626262"),
        legend.position = "none")
Carotenoids
## Saved with 4x3 inches as PDF