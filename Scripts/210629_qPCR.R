library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(Hmisc)
library(ggsci)
library(gridExtra)
library(readr)
library(stringr)
library(ggpubr)
library(extrafont)

Plates <- paste0("Plate", seq(1:5))

for (plate in Plates){
  file <- grep("Replicate", list.files(path = paste0("G:/Mi unidad/CRAG/qPCR/210629/", plate, "/")), value = TRUE)
  assign(paste("qPCR", plate, sep = "_"), read.csv(paste0("G:/Mi unidad/CRAG/qPCR/210629/", plate, "/", file), comment.char = "#") %>%
           dplyr::select(Sample, Target, Cq.Mean) %>%
           spread(Target, Cq.Mean) %>%
           filter(Sample!= "CAL", Sample!="crtB_0h_1", Sample!="Water")
         )
}

qPCR_Plate4 <- qPCR_Plate4 %>% dplyr::select(-ACT)
qPCR_Efficiency <- qPCR_Plate5[1:4, c(1,5:7)]
qPCR_Plate5 <- qPCR_Plate5[5:33, c(1,4)]

qPCR_All <- merge(qPCR_Plate1, qPCR_Plate2, by= "Sample")
qPCR_All <- merge(qPCR_All, qPCR_Plate3, by= "Sample")
qPCR_All <- merge(qPCR_All, qPCR_Plate4, by= "Sample", no.dups = TRUE)
qPCR_All <- merge(qPCR_All, qPCR_Plate5, by= "Sample", no.dups = TRUE)


qPCR_All <- qPCR_All %>% dplyr::mutate(ACT= coalesce(ACT.x, ACT.y),
                                       Hsp21= coalesce(Hsp21.x, Hsp21.y),
                                       Chap60B= coalesce(Chap60B.x, Chap60B.y)) %>%
  dplyr::select(Sample, ACT, FIB2B, Hsp21, Hsp3, Hsp18, Chap60B) %>%
  separate(Sample, into= c("Construct", "TimePoint", "Replicate"), sep="_")

qPCR_dCT <- dplyr::mutate(qPCR_All,
                dCT_FIB2B = 2^ACT / 2^FIB2B,
                dCT_Hsp21 = 2^ACT / 2^Hsp21,
                dCT_Hsp3 = 2^ACT / 2^Hsp3,
                dCT_Hsp18 = 2^ACT / 2^Hsp18,
                dCT_Chap60B = 2^ACT / 2^Chap60B) %>%
  dplyr::select(Construct, TimePoint, Replicate, dCT_FIB2B, dCT_Hsp21, dCT_Hsp3, dCT_Hsp18, dCT_Chap60B) %>%
  gather(key = "Gene", value= "dCT", -Construct, -TimePoint, -Replicate)

qPCR_dCT$Gene <- str_sub(qPCR_dCT$Gene, 5, -1)
qPCR_dCT$Construct <- factor(qPCR_dCT$Construct, levels = c("GFP", "crtB"), labels = c("GFP", "pcrtB"))
qPCR_dCT$TimePoint <- factor(qPCR_dCT$TimePoint, levels = unique(qPCR_dCT$TimePoint))
qPCR_dCT$Replicate <- as.factor(qPCR_dCT$Replicate)

qPCR_plot <- ggplot(qPCR_dCT %>% 
                      filter(Gene %in% c("Chap60B", "Hsp18", "Hsp3")) %>%
                      filter(TimePoint != "0h"), aes(x= TimePoint, y= dCT, group= Construct))+
  geom_point(shape = 20, size = 3, aes(col = Construct), alpha = 0.5)+
  geom_line(data = qPCR_dCT %>% 
              filter(Gene %in% c("Chap60B", "Hsp18", "Hsp3")) %>%
              filter(TimePoint != "0h") %>%
              select(Gene, Construct, TimePoint, dCT) %>%
              group_by(Gene, Construct, TimePoint)%>%
              summarise(AVG = mean(dCT)),
            aes(x = TimePoint, y = AVG, col = Construct),
            size = 1)+
  geom_crossbar(data = qPCR_dCT %>% 
                  filter(Gene %in% c("Chap60B", "Hsp18", "Hsp3")) %>%
                  filter(TimePoint != "0h") %>%
                  select(Gene, Construct, TimePoint, dCT) %>%
                  group_by(Gene, Construct, TimePoint)%>%
                  summarise(AVG = mean(dCT)),
                aes(x = TimePoint, y = AVG, ymin = AVG, ymax = AVG, col = Construct))+
  scale_fill_manual(values = c("GFP"= "#66c2a5",
                               "pcrtB" = "#fc8d62"),
                    label= c("GFP", "(p)crtB"))+
  scale_color_manual(values = c("GFP"= "#1b9e77",
                                "pcrtB" = "#d95f02"),
                     label= c("GFP", "(p)crtB"))+
  scale_y_continuous("Relative expression to ACT")+
  scale_x_discrete(label= c("24 hpi", "48 hpi", "72 hpi", "96 hpi"))+
  labs(title = "Expression qPCR", x = "Time Points")+
  facet_wrap(~Gene, scales = "free_y", ncol = 1)+
  theme_HPLC+
  theme(axis.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 8),
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey"))

qPCR_plot

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





Niben_Primers_list <- list("FIB2"="Niben101Scf02953g00007.1",
                           "Hsp3"="Niben101Scf09262g01017.1",
                           "Hsp21"="Niben101Scf08206g01019.1",
                           "Chap60B"="Niben101Scf10400g00005.1",
                           "Hsp18"="Niben101Scf01475g00019.1"
                           )

Niben_Primers_df <- as.data.frame(Niben_Primers_list)
Niben_Primers_df <- Niben_Primers_df %>% gather(key= "Protein", value = "Name")

Niben_Primers_261 <- read.table("G:/Mi unidad/CRAG/RNAseq/Milan/RNAseq_Milan/Funcional Groups/Niben_Primers_261.txt")
Niben_Primers_261 <- Niben_Primers_261[order(Niben_Primers_261$V1,-Niben_Primers_261$V3),]
Niben_Primers_261 <- Niben_Primers_261[!duplicated(Niben_Primers_261$V1),] %>%
  dplyr::select(1,2)

Niben_Primers_261_2 <- merge(Niben_Primers_df, Niben_Primers_261, by.x= "Name", by.y = 1) %>%
  dplyr::select(3,2)
colnames(Niben_Primers_261_2) <- c("Name", "Protein")

TPM_All_FC_Annotated <- read_delim("G:/Mi unidad/CRAG/RNAseq/Milan/RNAseq_Milan/Summary_AllSamples_TPM_FC.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
TPM_All_FC_Annotated <- as.data.frame(TPM_All_FC_Annotated)

TPM_All_Samples <- TPM_All_FC_Annotated %>% dplyr::select(1, 3:50)
TPM_All_Samples <- TPM_All_Samples %>% dplyr::mutate(GFP_22h= GFP22h1+GFP22h2+GFP22h3) %>%
  dplyr::mutate(GFP_25h= GFP25h1+GFP25h2+GFP25h3) %>%
  dplyr::mutate(GFP_28h= GFP28h1+GFP28h2+GFP28h3) %>%
  dplyr::mutate(GFP_34h= GFP34h1+GFP34h2+GFP34h3) %>%
  dplyr::mutate(GFP_37h= GFP37h1+GFP37h2+GFP37h3) %>%
  dplyr::mutate(GFP_40h= GFP40h1+GFP40h2+GFP40h3) %>%
  dplyr::mutate(GFP_46h= GFP46h1+GFP46h2+GFP46h3) %>%
  dplyr::mutate(GFP_56h= GFP56h1+GFP56h2+GFP56h3) %>%
  dplyr::mutate(pcrtB_22h= pcrtB22h1+pcrtB22h2+pcrtB22h3) %>%
  dplyr::mutate(pcrtB_25h= pcrtB25h1+pcrtB25h2+pcrtB25h3) %>%
  dplyr::mutate(pcrtB_28h= pcrtB28h1+pcrtB28h2+pcrtB28h3) %>%
  dplyr::mutate(pcrtB_34h= pcrtB34h1+pcrtB34h2+pcrtB34h3) %>%
  dplyr::mutate(pcrtB_37h= pcrtB37h1+pcrtB37h2+pcrtB37h3) %>%
  dplyr::mutate(pcrtB_40h= pcrtB40h1+pcrtB40h2+pcrtB40h3) %>%
  dplyr::mutate(pcrtB_46h= pcrtB46h1+pcrtB46h2+pcrtB46h3) %>%
  dplyr::mutate(pcrtB_56h= pcrtB56h1+pcrtB56h2+pcrtB56h3)

TPM_All_Samples <- TPM_All_Samples %>% dplyr::mutate(FC_22h= log2(pcrtB_22h / GFP_22h)) %>%
  dplyr::mutate(FC_25h= log2(pcrtB_25h / GFP_25h)) %>%
  dplyr::mutate(FC_28h= log2(pcrtB_28h / GFP_28h)) %>%
  dplyr::mutate(FC_34h= log2(pcrtB_34h / GFP_34h)) %>%
  dplyr::mutate(FC_37h= log2(pcrtB_37h / GFP_37h)) %>%
  dplyr::mutate(FC_40h= log2(pcrtB_40h / GFP_40h)) %>%
  dplyr::mutate(FC_46h= log2(pcrtB_46h / GFP_46h)) %>%
  dplyr::mutate(FC_56h= log2(pcrtB_56h / GFP_56h))

Cornell_FC_All <- TPM_All_Samples %>% dplyr::select(1, 66:73)


## Representing Heatmap

Niben_Primers_Ann <- Niben_Primers_261_2 %>% select(Name, Protein)

## We merge it with the big df created in the "Summary_RNAseq_Milan.R script with all the information

Niben_Primers <- merge(Niben_Primers_Ann, TPM_All_FC_Annotated, by="Name")

## We select only the FC data from the different comparisons
Niben_Primers_TPM <- Niben_Primers %>% select(Protein, 4:51, 61:66)
Niben_Primers_TPM <- Niben_Primers_TPM[!duplicated(Niben_Primers_TPM$Protein),]
Niben_Primers_TPM <- Niben_Primers_TPM[order(match(Niben_Primers_TPM[,1], Niben_Primers_Ann[,2])),]


## We select only the FC data from the different comparisons
Niben_Primers_FC <- Niben_Primers %>% select(Protein, 71, 73:87)
Niben_Primers_FC <- Niben_Primers_FC[!duplicated(Niben_Primers_FC$Protein),]
Niben_Primers_FC <- Niben_Primers_FC[order(match(Niben_Primers_FC[,1], Niben_Primers_Ann[,2])),]


## We select only the FC data from the different comparisons
Niben_Primers_DEG <- Niben_Primers %>% select(Protein, 89, 91:105)
Niben_Primers_DEG <- Niben_Primers_DEG[!duplicated(Niben_Primers_DEG$Protein),]
Niben_Primers_DEG <- Niben_Primers_DEG[order(match(Niben_Primers_DEG[,1], Niben_Primers_Ann[,2])),]

Niben_Primers_DEG$pcrtB22h_vs_GFP22h_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtB22h_vs_GFP22h_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtB25h_vs_GFP25h_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtB25h_vs_GFP25h_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtBAA_vs_GFPAA_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtBAA_vs_GFPAA_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtBBB_vs_GFPBB_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtBBB_vs_GFPBB_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtBYY_vs_GFPYY_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtBYY_vs_GFPYY_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtB28h_vs_GFP28h_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtB28h_vs_GFP28h_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtB34h_vs_GFP34h_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtB34h_vs_GFP34h_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtB37h_vs_GFP37h_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtB37h_vs_GFP37h_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtBZZ_vs_GFPZZ_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtBZZ_vs_GFPZZ_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtBCC_vs_GFPCC_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtBCC_vs_GFPCC_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtB40h_vs_GFP40h_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtB40h_vs_GFP40h_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtBDD_vs_GFPDD_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtBDD_vs_GFPDD_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtB46h_vs_GFP46h_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtB46h_vs_GFP46h_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtB56h_vs_GFP56h_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtB56h_vs_GFP56h_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtBEE_vs_GFPEE_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtBEE_vs_GFPEE_filtered) == TRUE, " ", "*")
Niben_Primers_DEG$pcrtB96h_vs_GFP96h_filtered <- ifelse(is.na(Niben_Primers_DEG$pcrtB96h_vs_GFP96h_filtered) == TRUE, " ", "*")
row.names(Niben_Primers_DEG) <- Niben_Primers_DEG$Protein
Niben_Primers_DEG[1] <- NULL

Niben_Primers_DEG <- Niben_Primers_DEG %>%
  select(pcrtB22h_vs_GFP22h_filtered, pcrtBAA_vs_GFPAA_filtered, pcrtB25h_vs_GFP25h_filtered, 
         pcrtBYY_vs_GFPYY_filtered, pcrtBBB_vs_GFPBB_filtered, pcrtB28h_vs_GFP28h_filtered, pcrtB34h_vs_GFP34h_filtered, pcrtBZZ_vs_GFPZZ_filtered, 
         pcrtB37h_vs_GFP37h_filtered, pcrtBCC_vs_GFPCC_filtered, pcrtB40h_vs_GFP40h_filtered, pcrtBDD_vs_GFPDD_filtered, pcrtB46h_vs_GFP46h_filtered, 
         pcrtBEE_vs_GFPEE_filtered, pcrtB56h_vs_GFP56h_filtered, pcrtB96h_vs_GFP96h_filtered)

row.names(Niben_Primers_FC) <- Niben_Primers_FC$Protein
Niben_Primers_FC[1] <- NULL

Niben_Primers_FC <- Niben_Primers_FC %>%
  select(pcrtB22h_vs_GFP22h_All, pcrtBAA_vs_GFPAA_All, pcrtB25h_vs_GFP25h_All, 
         pcrtBYY_vs_GFPYY_All, pcrtBBB_vs_GFPBB_All, pcrtB28h_vs_GFP28h_All, pcrtB34h_vs_GFP34h_All, pcrtBZZ_vs_GFPZZ_All, 
         pcrtB37h_vs_GFP37h_All, pcrtBCC_vs_GFPCC_All, pcrtB40h_vs_GFP40h_All, pcrtBDD_vs_GFPDD_All, pcrtB46h_vs_GFP46h_All, 
         pcrtBEE_vs_GFPEE_All, pcrtB56h_vs_GFP56h_All, pcrtB96h_vs_GFP96h_All)


Niben_Primers_TPM_ggplot <- Niben_Primers_TPM
Niben_Primers_TPM_ggplot <- Niben_Primers_TPM_ggplot %>% gather(key= "Sample", value= "value", -Protein)
Niben_Primers_TPM_ggplot$Construct <- str_sub(Niben_Primers_TPM_ggplot$Sample, 1, -5)
Niben_Primers_TPM_ggplot$Construct[Niben_Primers_TPM_ggplot$Construct=="S19_GFP4d_1"] <- "GFP"
Niben_Primers_TPM_ggplot$Construct[Niben_Primers_TPM_ggplot$Construct=="S20_GFP4d_2"] <- "GFP"
Niben_Primers_TPM_ggplot$Construct[Niben_Primers_TPM_ggplot$Construct=="S21_GFP4d_4"] <- "GFP"
Niben_Primers_TPM_ggplot$Construct[Niben_Primers_TPM_ggplot$Construct=="S22_pcrtB4d_2"] <- "pcrtB"
Niben_Primers_TPM_ggplot$Construct[Niben_Primers_TPM_ggplot$Construct=="S23_pcrtB4d_3"] <- "pcrtB"
Niben_Primers_TPM_ggplot$Construct[Niben_Primers_TPM_ggplot$Construct=="S24_pcrtB4d_4"] <- "pcrtB"
Niben_Primers_TPM_ggplot$Construct <- as.factor(Niben_Primers_TPM_ggplot$Construct)
Niben_Primers_TPM_ggplot$TimePoint <- str_sub(Niben_Primers_TPM_ggplot$Sample, -4, -2)
Niben_Primers_TPM_ggplot$TimePoint[Niben_Primers_TPM_ggplot$TimePoint=="_S1"] <- "96h"
Niben_Primers_TPM_ggplot$TimePoint[Niben_Primers_TPM_ggplot$TimePoint=="_S2"] <- "96h"
Niben_Primers_TPM_ggplot$TimePoint <- as.factor(Niben_Primers_TPM_ggplot$TimePoint)
Niben_Primers_TPM_ggplot$Replicate <- str_sub(Niben_Primers_TPM_ggplot$Sample, -1)
Niben_Primers_TPM_ggplot$Replicate <- as.factor(Niben_Primers_TPM_ggplot$Replicate)

'%ni%' <- Negate('%in%')

Niben_Primers_TPM_ggplot$TimePoint <- factor(Niben_Primers_TPM_ggplot$TimePoint, 
                                             levels = c("25h", "46h", "72h", "96h"))

RNAseq_plot <- ggplot(Niben_Primers_TPM_ggplot %>% 
                        filter(Protein %in% c("Chap60B", "Hsp18", "Hsp3")) %>%
                        filter(TimePoint %in% c("25h", "46h", "72h", "96h")), aes(x = TimePoint, y = value, group = Construct))+
  geom_point(shape = 20, size = 3, aes(col = Construct), alpha = 0.5)+
  #stat_summary(fun.data = mean_sdl, geom = "errorbar",fun.args = list(mult = 1), width = 0.4)+
  geom_line(data = Niben_Primers_TPM_ggplot %>% 
              filter(Protein %in% c("Chap60B", "Hsp18", "Hsp3")) %>%
              filter(TimePoint %in% c("25h", "46h", "72h", "96h")) %>%
              select(Protein, Construct, TimePoint, value) %>%
              group_by(Protein, Construct, TimePoint)%>%
              summarise(AVG = mean(value)),
            aes(x = TimePoint, y = AVG, col = Construct),
            size = 1)+
  geom_crossbar(data = Niben_Primers_TPM_ggplot %>% 
                  filter(Protein %in% c("Chap60B", "Hsp18", "Hsp3")) %>%
                  filter(TimePoint %in% c("25h", "46h", "72h", "96h")) %>%
                  select(Protein, Construct, TimePoint, value) %>%
                  group_by(Protein, Construct, TimePoint)%>%
                  summarise(AVG = mean(value)),
                aes(x = TimePoint, y = AVG, ymin = AVG, ymax = AVG, col = Construct))+
  scale_fill_manual(values = c("GFP"= "#66c2a5",
                               "pcrtB" = "#fc8d62"),
                    label= c("GFP", "(p)crtB"))+
  scale_color_manual(values = c("GFP"= "#1b9e77",
                                "pcrtB" = "#d95f02"),
                     label= c("GFP", "(p)crtB"))+
  scale_y_continuous("TPM values")+
  scale_x_discrete(limits= c("25h", "46h", "72h", "96h"),
                   label = c("25 hpi", "46 hpi", "72 hpi", "96 hpi"))+
  labs(title = "Expression RNAseq", x = "Time Points")+
  facet_wrap(~Protein, scales = "free_y", ncol = 1)+
  theme_HPLC+
  theme(axis.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        title = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 8),
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        panel.grid.major.y = element_line(linetype = "dotted", colour = "grey"))

ggarrange(RNAseq_plot, qPCR_plot, ncol = 2, nrow = 1)


theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
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
