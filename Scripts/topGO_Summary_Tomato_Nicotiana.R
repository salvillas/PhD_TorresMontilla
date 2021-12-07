library(DESeq2)
library(ggplot2)
library(gplots)
library(genefilter)
library(GenomicRanges)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(Biobase)
library(gplots)
library(graphics)
library(grDevices)
library(survival)
library(UpSetR)
library(ComplexHeatmap)
library(topGO)
library(ggsci)
library(scales)
library(tidyverse)
library(ggpubr)

###################################################################################################
## This is a Script to summarize all GO terms from Nicotiana & Tomato & Arabidopsis experiments
###################################################################################################

##### First we work with the Shinozaki experiment

setwd("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Tomato/2nd/Shinozaki/06_GOterms")

## We create a list with all files of the folder from REVIGO
sampleFiles_Shinozaki <- grep("REVIGO", list.files("."), value = TRUE)

## We load all files with a for loop, adding a column with the experiment label
## and another with the background. We select only those with 0.7 dispensability or less
for (element in sampleFiles_Shinozaki){
  object1 <- str_remove(element, ".csv")
  assign(object1, read.csv(element, header = TRUE))
  nam <- assign(object1, read.csv(element, header = TRUE))
  Set <- "Shinozaki"
  nam$Set <- Set
  Group <- str_remove(object1, "topGO_weightFS_")
  Group <- str_remove(Group, "_BP_4.0_AHRD_REVIGO")
  nam$Group <- Group
  Background <- "M82"
  nam$Background <- Background
  nam <- nam %>% 
    dplyr::filter(dispensability<=0.7) %>%
    dplyr::select(description, Group, Set, Background) %>%
    tidyr::separate(Group, into=c("Stage", "UpDown"), sep="_")
  assign(object1, nam)
}

##### Now the same with Yazdani experiment

setwd("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Tomato/2nd/Yazdani/GO_Terms")

sampleFiles_Yazdani <- grep("REVIGO", list.files("."), value = TRUE)

## We repeat the for loop, but now the background comes from the file name
for (element in sampleFiles_Yazdani){
  object1 <- str_remove(element, ".csv")
  assign(object1, read.csv(element, header = TRUE))
  nam <- assign(object1, read.csv(element, header = TRUE))
  Set <- "Yazdani"
  nam$Set <- Set
  Stage <- str_remove(object1, "topGO_weightFS_")
  Stage <- str_remove(Stage, "_BP_4.0_AHRD_REVIGO")
  nam$Stage <- Stage
  nam <- nam %>% 
    dplyr::filter(dispensability<=0.7) %>%
    dplyr::select(description, Stage, Set) %>%
    tidyr::separate(Stage, into=c("Background", "Stage", "UpDown"), sep="_") %>%
    dplyr::select(description, Stage, UpDown, Set, Background)
  assign(object1, nam)
}


##### Now the same with Nicotiana experiments

setwd("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd")

sampleFiles_Niben <- grep("REVIGO.csv", list.files("."), value = TRUE)

## We repeat the for loop, but now the background comes from the file name
for (element in sampleFiles_Niben){
  object1 <- str_remove(element, ".csv")
  assign(object1, read.csv(element, header = TRUE))
  nam <- assign(object1, read.csv(element, header = TRUE))
  Set <- "Nicotiana"
  nam$Set <- Set
  Stage <- str_remove(object1, "topGO_weightFS_")
  Stage <- str_remove(Stage, "_BP_REVIGO")
  nam$Stage <- Stage
  Background <- "Niben"
  nam$Background <- Background
  nam <- nam %>% 
    dplyr::filter(dispensability<=0.7) %>%
    dplyr::select(description, Stage, Set, Background) %>%
    tidyr::separate(Stage, into=c("Stage", "UpDown"), sep="_") %>%
    dplyr::select(description, Stage, UpDown, Set, Background)
  assign(object1, nam)
}

##### Now the same with Pepper experiment

setwd("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Pepper")

sampleFiles_Pepper <- grep("REVIGO.csv", list.files("."), value = TRUE)

## We repeat the for loop, but now the background comes from the file name
for (element in sampleFiles_Pepper){
  object1 <- str_remove(element, ".csv")
  assign(object1, read.csv(element, header = TRUE))
  nam <- assign(object1, read.csv(element, header = TRUE))
  Set <- "Capsicum"
  nam$Set <- Set
  Stage <- str_remove(object1, "topGO_weightFS_")
  Stage <- str_remove(Stage, "_BP_REVIGO")
  nam$Stage <- Stage
  Background <- "Pepper"
  nam$Background <- Background
  nam <- nam %>% 
    dplyr::filter(dispensability<=0.7) %>%
    dplyr::select(description, Stage, Set, Background) %>%
    tidyr::separate(Stage, into=c("Stage", "UpDown"), sep="_") %>%
    dplyr::select(description, Stage, UpDown, Set, Background)
  assign(object1, nam)
}

##### Now the same with Arabidopsis experiment

setwd("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Arabidopsis/")

sampleFiles_Arabidopsis <- grep("REVIGO.csv", list.files("."), value = TRUE)

## We repeat the for loop, but now the background comes from the file name
for (element in sampleFiles_Arabidopsis){
  object1 <- str_remove(element, ".csv")
  assign(object1, read.csv(element, header = TRUE))
  nam <- assign(object1, read.csv(element, header = TRUE))
  Set <- "Arabidopsis"
  nam$Set <- Set
  Stage <- str_remove(object1, "topGO_weightFS_")
  Stage <- str_remove(Stage, "_BP_REVIGO")
  nam$Stage <- Stage
  Background <- "Ath"
  nam$Background <- Background
  nam <- nam %>% 
    dplyr::filter(Dispensability<=0.7) %>%
    dplyr::select(Name, Stage, Set, Background) %>%
    tidyr::separate(Stage, into=c("Stage", "UpDown"), sep="_") %>%
    dplyr::select(Name, Stage, UpDown, Set, Background)
  nam$Name <- str_sub(nam$Name, 2, -1)
  colnames(nam) <- c("description", "Stage", "UpDown", "Set", "Background")
  assign(object1, nam)
}

##################################################################################################
########################################  TIME POINTS  ###########################################
##################################################################################################
## YAZDANI VS NICOTIANA

########### Up

filelist_Yazdani_vs_Niben_Up <- lapply(list("topGO_weightFS_M82_B_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_M82_or_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_M82_Red_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_MG_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_B_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_22h_Up_BP_REVIGO",
                                            "topGO_weightFS_25h_Up_BP_REVIGO",
                                            "topGO_weightFS_28h_Up_BP_REVIGO",
                                            "topGO_weightFS_34h_Up_BP_REVIGO",
                                            "topGO_weightFS_37h_Up_BP_REVIGO",
                                            "topGO_weightFS_40h_Up_BP_REVIGO",
                                            "topGO_weightFS_46h_Up_BP_REVIGO",
                                            "topGO_weightFS_56h_Up_BP_REVIGO",
                                            "topGO_weightFS_96h_Up_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Niben_Up <- do.call(rbind, filelist_Yazdani_vs_Niben_Up)

Summary_GOterms_Yazdani_vs_Niben_Up$Stage <- factor(Summary_GOterms_Yazdani_vs_Niben_Up$Stage,
                                                    levels = c("B", "or", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h"),
                                                    labels = c("Breaker", "Orange", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Niben_Up$Background <- factor(Summary_GOterms_Yazdani_vs_Niben_Up$Background,
                                                         levels = c("M82", "MG", "B", "Niben"),
                                                         labels = c("M82", "MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Niben_Up$Group <- paste(Summary_GOterms_Yazdani_vs_Niben_Up$Stage, Summary_GOterms_Yazdani_vs_Niben_Up$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Niben_Up$Group <- factor(Summary_GOterms_Yazdani_vs_Niben_Up$Group,
                                                    levels = c("Breaker_M82", "Orange_M82", "Red_M82", "ORHis_MG", "ORHis_Breaker","22h_Niben", "25h_Niben", "28h_Niben", "34h_Niben", "37h_Niben", "40h_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Niben_Up <- Summary_GOterms_Yazdani_vs_Niben_Up[order(Summary_GOterms_Yazdani_vs_Niben_Up$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_Up)){
  Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Niben_Up$Group[Summary_GOterms_Yazdani_vs_Niben_Up$description==Summary_GOterms_Yazdani_vs_Niben_Up$description[i]])
  Summary_GOterms_Yazdani_vs_Niben_Up$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Niben_Up$Group[Summary_GOterms_Yazdani_vs_Niben_Up$description==Summary_GOterms_Yazdani_vs_Niben_Up$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_Up)){
  if(Summary_GOterms_Yazdani_vs_Niben_Up$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Niben_Up$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Niben_Up$color[i] <- "#8DA0CB"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Up$NumberSets[i]>11){
    Summary_GOterms_Yazdani_vs_Niben_Up$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Niben_Up$color[i] <- "#66C2A5"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Niben_Up$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Niben_Up$color[i] <- "#E78AC3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Niben_Up$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Niben_Up$color[i] <- "#FC8D62"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i]))>5)){
    Summary_GOterms_Yazdani_vs_Niben_Up$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Niben_Up$color[i] <- "#B3B3B3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i]))>3) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i])) %in% c(6:14))){
    Summary_GOterms_Yazdani_vs_Niben_Up$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Niben_Up$color[i] <- "#FFD92F"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Up$TotalSets[i])) %in% c(1:3))){
    Summary_GOterms_Yazdani_vs_Niben_Up$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Niben_Up$color[i] <- "#E5C494"
  } else {
    Summary_GOterms_Yazdani_vs_Niben_Up$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Niben_Up$color[i] <- "#A6D854"}
}

## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Niben_Up$freq <- factor(Summary_GOterms_Yazdani_vs_Niben_Up$freq, levels = c("Constant", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Niben_Up$description <- factor(Summary_GOterms_Yazdani_vs_Niben_Up$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Niben_Up$description[order(Summary_GOterms_Yazdani_vs_Niben_Up$Group, Summary_GOterms_Yazdani_vs_Niben_Up$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors1 <- Summary_GOterms_Yazdani_vs_Niben_Up %>% dplyr::select(description, color)
colors1 <- unique(colors1)
colors1 <- colors1[order(colors1$description),]
colors1 <- colors1$color
names(colors1) <- rev(unique(Summary_GOterms_Yazdani_vs_Niben_Up$description[order(Summary_GOterms_Yazdani_vs_Niben_Up$Group, Summary_GOterms_Yazdani_vs_Niben_Up$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Niben_Up_plot <- ggplot(Summary_GOterms_Yazdani_vs_Niben_Up, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#66C2A5",
                               "Ripening"="#E78AC3",
                               "Nicotiana-independent"="#E5C494",
                               "OR-Dependent"="#FC8D62", 
                               "Ripening-independent"="#FFD92F", 
                               "Spread"="#A6D854", 
                               "Nicotiana"="#B3B3B3",
                               "Exclusive"="#8DA0CB"),
                     breaks = c("Constant",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Nicotiana Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme(axis.text.y = element_text(color = colors1, face = "bold"))
Summary_GOterms_Yazdani_vs_Niben_Up_plot ## Saved pdf 20x70 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Yazdani_vs_Niben_Up %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Yazdani_vs_Niben_Up_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


########### Down

filelist_Yazdani_vs_Niben_Down <- lapply(list("topGO_weightFS_M82_B_Down_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_M82_or_Down_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_M82_Red_Down_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_MG_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_B_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_22h_Down_BP_REVIGO",
                                            "topGO_weightFS_25h_Down_BP_REVIGO",
                                            "topGO_weightFS_28h_Down_BP_REVIGO",
                                            "topGO_weightFS_34h_Down_BP_REVIGO",
                                            "topGO_weightFS_37h_Down_BP_REVIGO",
                                            "topGO_weightFS_40h_Down_BP_REVIGO",
                                            "topGO_weightFS_46h_Down_BP_REVIGO",
                                            "topGO_weightFS_56h_Down_BP_REVIGO",
                                            "topGO_weightFS_96h_Down_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Niben_Down <- do.call(rbind, filelist_Yazdani_vs_Niben_Down)

Summary_GOterms_Yazdani_vs_Niben_Down$Stage <- factor(Summary_GOterms_Yazdani_vs_Niben_Down$Stage,
                                                    levels = c("B", "or", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h"),
                                                    labels = c("Breaker", "Orange", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Niben_Down$Background <- factor(Summary_GOterms_Yazdani_vs_Niben_Down$Background,
                                                         levels = c("M82", "MG", "B", "Niben"),
                                                         labels = c("M82", "MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Niben_Down$Group <- paste(Summary_GOterms_Yazdani_vs_Niben_Down$Stage, Summary_GOterms_Yazdani_vs_Niben_Down$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Niben_Down$Group <- factor(Summary_GOterms_Yazdani_vs_Niben_Down$Group,
                                                    levels = c("Breaker_M82", "Orange_M82", "Red_M82", "ORHis_MG", "ORHis_Breaker","22h_Niben", "25h_Niben", "28h_Niben", "34h_Niben", "37h_Niben", "40h_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Niben_Down <- Summary_GOterms_Yazdani_vs_Niben_Down[order(Summary_GOterms_Yazdani_vs_Niben_Down$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_Down)){
  Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Niben_Down$Group[Summary_GOterms_Yazdani_vs_Niben_Down$description==Summary_GOterms_Yazdani_vs_Niben_Down$description[i]])
  Summary_GOterms_Yazdani_vs_Niben_Down$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Niben_Down$Group[Summary_GOterms_Yazdani_vs_Niben_Down$description==Summary_GOterms_Yazdani_vs_Niben_Down$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_Down)){
  if(Summary_GOterms_Yazdani_vs_Niben_Down$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Niben_Down$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Niben_Down$color[i] <- "#8DA0CB"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Down$NumberSets[i]>11){
    Summary_GOterms_Yazdani_vs_Niben_Down$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Niben_Down$color[i] <- "#66C2A5"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Niben_Down$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Niben_Down$color[i] <- "#E78AC3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Niben_Down$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Niben_Down$color[i] <- "#FC8D62"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i]))>5)){
    Summary_GOterms_Yazdani_vs_Niben_Down$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Niben_Down$color[i] <- "#B3B3B3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i]))>3) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i])) %in% c(6:14))){
    Summary_GOterms_Yazdani_vs_Niben_Down$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Niben_Down$color[i] <- "#FFD92F"
  } else if(Summary_GOterms_Yazdani_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_Down$TotalSets[i])) %in% c(1:3))){
    Summary_GOterms_Yazdani_vs_Niben_Down$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Niben_Down$color[i] <- "#E5C494"
  } else {
    Summary_GOterms_Yazdani_vs_Niben_Down$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Niben_Down$color[i] <- "#A6D854"}
}

## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Niben_Down$freq <- factor(Summary_GOterms_Yazdani_vs_Niben_Down$freq, levels = c("Constant", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Niben_Down$description <- factor(Summary_GOterms_Yazdani_vs_Niben_Down$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Niben_Down$description[order(Summary_GOterms_Yazdani_vs_Niben_Down$Group, Summary_GOterms_Yazdani_vs_Niben_Down$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors2 <- Summary_GOterms_Yazdani_vs_Niben_Down %>% dplyr::select(description, color)
colors2 <- unique(colors2)
colors2 <- colors2[order(colors2$description),]
colors2 <- colors2$color
names(colors2) <- rev(unique(Summary_GOterms_Yazdani_vs_Niben_Down$description[order(Summary_GOterms_Yazdani_vs_Niben_Down$Group, Summary_GOterms_Yazdani_vs_Niben_Down$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Niben_Down_plot <- ggplot(Summary_GOterms_Yazdani_vs_Niben_Down, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#66C2A5",
                               "Ripening"="#E78AC3",
                               "Nicotiana-independent"="#E5C494",
                               "OR-Dependent"="#FC8D62", 
                               "Ripening-independent"="#FFD92F", 
                               "Spread"="#A6D854", 
                               "Nicotiana"="#B3B3B3",
                               "Exclusive"="#8DA0CB"),
                     breaks = c("Constant",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Nicotiana Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme(axis.text.y = element_text(color = colors2, face = "bold"))
Summary_GOterms_Yazdani_vs_Niben_Down_plot ## Saved pdf 20x70 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Yazdani_vs_Niben_Down %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Yazdani_vs_Niben_Down_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)



##################################################################################################
## YAZDANI VS PEPPER VS NICOTIANA

########### Up

filelist_Yazdani_vs_Pepper_vs_Niben_Up <- lapply(list("topGO_weightFS_M82_B_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_M82_or_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_M82_Red_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_Breaker_Up_BP_REVIGO",
                                                      "topGO_weightFS_Red_Up_BP_REVIGO",
                                                      "topGO_weightFS_MG_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_B_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_22h_Up_BP_REVIGO",
                                                      "topGO_weightFS_25h_Up_BP_REVIGO",
                                                      "topGO_weightFS_28h_Up_BP_REVIGO",
                                                      "topGO_weightFS_34h_Up_BP_REVIGO",
                                                      "topGO_weightFS_37h_Up_BP_REVIGO",
                                                      "topGO_weightFS_40h_Up_BP_REVIGO",
                                                      "topGO_weightFS_46h_Up_BP_REVIGO",
                                                      "topGO_weightFS_56h_Up_BP_REVIGO",
                                                      "topGO_weightFS_96h_Up_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up <- do.call(rbind, filelist_Yazdani_vs_Pepper_vs_Niben_Up)

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Stage <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Stage,
                                                              levels = c("B", "or", "Red", "Breaker", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h"),
                                                              labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Background <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Background,
                                                                   levels = c("M82", "Pepper","MG", "B", "Niben"),
                                                                   labels = c("M82", "Pepper","MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Group <- paste(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Stage, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Group <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Group,
                                                              levels = c("Breaker_M82", "Orange_M82", "Red_M82", "Breaker_Pepper", "Red_Pepper","ORHis_MG", "ORHis_Breaker","22h_Niben", "25h_Niben", "28h_Niben", "34h_Niben", "37h_Niben", "40h_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up)){
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$description[i]])
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up)){
  if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i]>13){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#370335B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "Tomato-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#bea3ce"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "Pepper-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#fccde5"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i])) %in% c(1:3)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i])) %in% c(4,5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#bc80bd"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i]))>5) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i]))<8)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#F05C3BB2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i]))>7)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#D2AF81B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i]))>5) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i])) %in% c(8:16))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#FED439B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i]))<8) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$TotalSets[i])) %in% c(1:5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#808080"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq, levels = c("Constant","Tomato-Ripening", "Pepper-Ripening", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors3 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up %>% dplyr::select(description, color)
colors3 <- unique(colors3)
colors3 <- colors3[order(colors3$description),]
colors3 <- colors3$color
names(colors3) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#370335B2",
                               "Tomato-Ripening" = "#bea3ce",#"#bebada",
                               "Pepper-Ripening" = "#fccde5",
                               "Ripening"="#bc80bd",
                               "Nicotiana-independent"="#808080",#"#D2AF81B2",
                               "OR-Dependent"="#F05C3BB2", 
                               "Ripening-independent"="#FED439B2", 
                               "Spread"="#b3de69", 
                               "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                               "Exclusive"="#197EC0B2"),
                     breaks = c("Constant",
                                "Tomato-Ripening",
                                "Pepper-Ripening",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme(axis.text.y = element_text(color = colors3, face = "bold"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_plot ## Saved pdf 23x73 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

########### Down

filelist_Yazdani_vs_Pepper_vs_Niben_Down <- lapply(list("topGO_weightFS_M82_B_Down_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_M82_or_Down_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_M82_Red_Down_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_Breaker_Down_BP_REVIGO",
                                                      "topGO_weightFS_Red_Down_BP_REVIGO",
                                                      "topGO_weightFS_MG_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_B_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_22h_Down_BP_REVIGO",
                                                      "topGO_weightFS_25h_Down_BP_REVIGO",
                                                      "topGO_weightFS_28h_Down_BP_REVIGO",
                                                      "topGO_weightFS_34h_Down_BP_REVIGO",
                                                      "topGO_weightFS_37h_Down_BP_REVIGO",
                                                      "topGO_weightFS_40h_Down_BP_REVIGO",
                                                      "topGO_weightFS_46h_Down_BP_REVIGO",
                                                      "topGO_weightFS_56h_Down_BP_REVIGO",
                                                      "topGO_weightFS_96h_Down_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down <- do.call(rbind, filelist_Yazdani_vs_Pepper_vs_Niben_Down)

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Stage <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Stage,
                                                              levels = c("B", "or", "Red", "Breaker", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h"),
                                                              labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Background <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Background,
                                                                   levels = c("M82", "Pepper","MG", "B", "Niben"),
                                                                   labels = c("M82", "Pepper","MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Group <- paste(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Stage, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Group <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Group,
                                                              levels = c("Breaker_M82", "Orange_M82", "Red_M82", "Breaker_Pepper", "Red_Pepper","ORHis_MG", "ORHis_Breaker","22h_Niben", "25h_Niben", "28h_Niben", "34h_Niben", "37h_Niben", "40h_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down)){
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$description[i]])
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down)){
  if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i]>13){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#370335B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "Tomato-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#bea3ce"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "Pepper-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#fccde5"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i])) %in% c(1:3)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i])) %in% c(4,5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#bc80bd"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i]))>5) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i]))<8)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#F05C3BB2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i]))>7)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#D2AF81B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i]))>5) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i])) %in% c(8:16))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#FED439B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i]))<8) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$TotalSets[i])) %in% c(1:5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#808080"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq, levels = c("Constant","Tomato-Ripening", "Pepper-Ripening", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors4 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down %>% dplyr::select(description, color)
colors4 <- unique(colors4)
colors4 <- colors4[order(colors4$description),]
colors4 <- colors4$color
names(colors4) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#370335B2",
                               "Tomato-Ripening" = "#bea3ce",#"#bebada",
                               "Pepper-Ripening" = "#fccde5",
                               "Ripening"="#bc80bd",
                               "Nicotiana-independent"="#808080",#"#D2AF81B2",
                               "OR-Dependent"="#F05C3BB2", 
                               "Ripening-independent"="#FED439B2", 
                               "Spread"="#b3de69", 
                               "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                               "Exclusive"="#197EC0B2"),
                     breaks = c("Constant",
                                "Tomato-Ripening",
                                "Pepper-Ripening",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme(axis.text.y = element_text(color = colors4, face = "bold"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_plot ## Saved pdf 23x73 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

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

######################################################
## SUMMARIZING RESULTS BY GROUP (% AND COUNTS)

## Up
percentData_Yazdani_vs_Pepper_vs_Niben_Up <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up %>%
  group_by(Group) %>% count(freq) %>%
  mutate(ratio=scales::percent(n/sum(n)))
percentData_Yazdani_vs_Pepper_vs_Niben_Up

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_percent <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up, aes(x= Group))+
  geom_bar(aes(fill=freq), position = position_fill(reverse = TRUE))+
  geom_text(data = percentData_Yazdani_vs_Pepper_vs_Niben_Up, aes(y=n, label=ratio),
            position = position_fill(vjust = 0.5), colour = "white", fontface = "bold")+
  scale_fill_manual(values= c("Constant"="#370335B2",
                              "Tomato-Ripening" = "#bea3ce",#"#bebada",
                              "Pepper-Ripening" = "#fccde5",
                              "Ripening"="#bc80bd",
                              "Nicotiana-independent"="#808080",#"#D2AF81B2",
                              "OR-Dependent"="#F05C3BB2", 
                              "Ripening-independent"="#FED439B2", 
                              "Spread"="#b3de69", 
                              "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                              "Exclusive"="#197EC0B2"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     labels = c("0", "25", "50", "75", "100"))+
  ylab("% of GO terms subgroups")+
  theme_HPLC+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_percent 

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_abs <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up, aes(x=Group))+
  geom_bar(aes(fill=Background))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_discrete(labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h"))+
  ylab("Number of GO terms")+
  scale_fill_manual(values = c("M82" = "#bea3ce",
                               "Pepper" ="#fccde5",
                               "MG"="#e45132",
                               "Breaker"="#bd2e15",
                               "Niben"="#D2AF81B2"))+
  theme_HPLC+
  theme(legend.position = "top",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_abs

figureA <- ggarrange(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_abs,
                     Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_percent,
                     ncol = 1,
                     nrow = 2)
figureA ## Saved pdf 15x10 "Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Percent"


## Down
percentData_Yazdani_vs_Pepper_vs_Niben_Down <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down %>%
  group_by(Group) %>% count(freq) %>%
  mutate(ratio=scales::percent(n/sum(n)))
percentData_Yazdani_vs_Pepper_vs_Niben_Down

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_percent <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down, aes(x= Group))+
  geom_bar(aes(fill=freq), position = position_fill(reverse = TRUE))+
  geom_text(data = percentData_Yazdani_vs_Pepper_vs_Niben_Down, aes(y=n, label=ratio),
            position = position_fill(vjust = 0.5), colour = "white", fontface = "bold")+
  scale_fill_manual(values= c("Constant"="#370335B2",
                              "Tomato-Ripening" = "#bea3ce",#"#bebada",
                              "Pepper-Ripening" = "#fccde5",
                              "Ripening"="#bc80bd",
                              "Nicotiana-independent"="#808080",#"#D2AF81B2",
                              "OR-Dependent"="#F05C3BB2", 
                              "Ripening-independent"="#FED439B2", 
                              "Spread"="#b3de69", 
                              "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                              "Exclusive"="#197EC0B2"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     labels = c("0", "25", "50", "75", "100"))+
  ylab("% of GO terms subGroups")+
  theme_HPLC+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_percent 

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_abs <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down, aes(x=Group))+
  geom_bar(aes(fill=Background))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_discrete(labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h"))+
  scale_fill_manual(values = c("M82" = "#bea3ce",
                               "Pepper" ="#fccde5",
                               "MG"="#e45132",
                               "Breaker"="#bd2e15",
                               "Niben"="#D2AF81B2"))+
  theme_HPLC+
  theme(legend.position = "top",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_abs

figureB <- ggarrange(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_abs,
                     Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_percent,
                     ncol = 1,
                     nrow = 2)
figureB ## Saved pdf 15x10 "Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Percent"



############################################################################
## Graph only Spread

## Up

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq=="Spread",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$TotalSets[i]))<13) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$TotalSets[i])) %in% c(1:5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$TotalSets[i])) %in% c(8:12))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$TotalSets[i]))%in% c(1:5)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$TotalSets[i])) %ni% c(8:12)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$TotalSets[i])) %in% c(13:16))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$Group)]))))

colors5 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread %>% dplyr::select(description, color2)
colors5 <- unique(colors5)
colors5 <- colors5[order(colors5$description),]
colors5 <- colors5$color2
names(colors5) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 16.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,16,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Spread Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors5, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread_plot ## Saved pdf 23x73 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Spread %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbTimePoints_Up_0.7_Spread.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


## Down

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq=="Spread",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$TotalSets[i]))<13) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$TotalSets[i])) %in% c(1:5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$TotalSets[i])) %in% c(8:12))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$TotalSets[i]))%in% c(1:5)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$TotalSets[i])) %ni% c(8:12)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$TotalSets[i])) %in% c(13:16))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$Group)]))))

colors6 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread %>% dplyr::select(description, color2)
colors6 <- unique(colors6)
colors6 <- colors6[order(colors6$description),]
colors6 <- colors6$color2
names(colors6) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 16.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,16,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Spread Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors6, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread_plot ## Saved pdf 23x73 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Spread %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbTimePoints_Down_0.7_Spread.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


############################################################################
## Graph only Ripening-independet

## Up

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up$freq=="Ripening-independent",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$TotalSets[i]))<13)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$TotalSets[i]))%in% c(6,7)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$TotalSets[i])) %ni% c(8:12)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$TotalSets[i])) %in% c(13:16))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$Group)]))))

colors7 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind %>% dplyr::select(description, color2)
colors7 <- unique(colors7)
colors7 <- colors7[order(colors7$description),]
colors7 <- colors7$color2
names(colors7) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 16.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,16,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Ripening-Independent Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors7, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind_plot ## Saved pdf 23x73 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Rip_ind %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbTimePoints_Up_0.7_Ripening_independent.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Down

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down$freq=="Ripening-independent",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$TotalSets[i]))<13)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$TotalSets[i]))%in% c(6,7)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$TotalSets[i])) %ni% c(8:12)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$TotalSets[i])) %in% c(13:16))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$Group)]))))

colors8 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind %>% dplyr::select(description, color2)
colors8 <- unique(colors8)
colors8 <- colors8[order(colors8$description),]
colors8 <- colors8$color2
names(colors8) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 16.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,16,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Ripening-Independent Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors8, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind_plot ## Saved pdf 23x73 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Rip_ind %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbTimePoints_Down_0.7_Ripening_independent.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


##################################################################################################
## SHINOZAKI VS YAZDANI VS PEPPER VS NICOTIANA VS ARABIDOPSIS

########### Up

filelist_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up <- lapply(list("topGO_weightFS_Pk_Up_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_LR_Up_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_RR_Up_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_M82_B_Up_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_M82_or_Up_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_M82_Red_Up_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_Breaker_Up_BP_REVIGO",
                                                                          "topGO_weightFS_Red_Up_BP_REVIGO",
                                                                          "topGO_weightFS_MG_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_B_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_22h_Up_BP_REVIGO",
                                                                          "topGO_weightFS_25h_Up_BP_REVIGO",
                                                                          "topGO_weightFS_28h_Up_BP_REVIGO",
                                                                          "topGO_weightFS_34h_Up_BP_REVIGO",
                                                                          "topGO_weightFS_37h_Up_BP_REVIGO",
                                                                          "topGO_weightFS_40h_Up_BP_REVIGO",
                                                                          "topGO_weightFS_46h_Up_BP_REVIGO",
                                                                          "topGO_weightFS_56h_Up_BP_REVIGO",
                                                                          "topGO_weightFS_96h_Up_BP_REVIGO",
                                                                          "topGO_weightFS_Woo_Up_BP_REVIGO"), get)
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up <- do.call(rbind, filelist_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up)

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Stage <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Stage,
                                                              levels = c("Pk", "LR", "RR", "B", "or", "Red", "Breaker", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h", "Woo"),
                                                              labels = c("Breaker", "Orange", "Red", "Breaker", "Orange", "Red", "Breaker", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h", "Senescence"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Background <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Background,
                                                                   levels = c("M82", "Pepper","MG", "B", "Niben", "Ath"),
                                                                   labels = c("M82", "Pepper","MG", "Breaker", "Niben", "Ath"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group <- paste(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Stage, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Background, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Set, sep = "_")
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group,
                                                              levels = c("Breaker_M82_Shinozaki", "Orange_M82_Shinozaki", "Red_M82_Shinozaki",
                                                                         "Breaker_M82_Yazdani", "Orange_M82_Yazdani", "Red_M82_Yazdani",
                                                                         "Breaker_Pepper_Capsicum", "Red_Pepper_Capsicum","ORHis_MG_Yazdani", "ORHis_Breaker_Yazdani",
                                                                         "22h_Niben_Nicotiana", "25h_Niben_Nicotiana", "28h_Niben_Nicotiana", "34h_Niben_Nicotiana", 
                                                                         "37h_Niben_Nicotiana", "40h_Niben_Nicotiana", "46h_Niben_Nicotiana", "56h_Niben_Nicotiana", 
                                                                         "96h_Niben_Nicotiana", "Senescence_Ath_Arabidopsis"),
                                                              labels = c("Breaker_Shinozaki", "Orange_Shinozaki", "Red_Shinozaki",
                                                                         "Breaker_Yazdani", "Orange_Yazdani", "Red_Yazdani",
                                                                         "Breaker_Pepper", "Red_Pepper","ORHis_MG", "ORHis_Breaker",
                                                                         "22h_Niben", "25h_Niben", "28h_Niben", "34h_Niben", 
                                                                         "37h_Niben", "40h_Niben", "46h_Niben", "56h_Niben", 
                                                                         "96h_Niben", "Senescence_Ath"))

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group),]

for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up)){
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i] <- list(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group[Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description==Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description[i]])
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i] <- length(c(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group[Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description==Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description[i]]))
}

'%ni%' <- Negate('%in%')
# First approach
for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up)){
  if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]==1){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Exclusive"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>13){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Constant"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#370335B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))<7)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Tomato-Ripening"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#bea3ce"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))>6) & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))<9)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Pepper-Ripening"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#fccde5"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))<9) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(1:6)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(7,8))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Ripening"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#bc80bd"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))>8) & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))<11)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "OR-Dependent"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#F05C3BB2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))>10) & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))<20)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Nicotiana"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#D2AF81B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))>8) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(9,10)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))>10)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Ripening-independent"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#FED439B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(1:10)) & all(!as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(11:19))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#808080"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))>10) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(11:19)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i]))>19)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Fruit-independent"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#006c00"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Spread"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq, levels = c("Constant","Tomato-Ripening", "Pepper-Ripening", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Fruit-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description, levels = (rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors3 <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up %>% dplyr::select(description, color)
colors3 <- unique(colors3)
colors3 <- colors3[order(colors3$description),]
colors3 <- colors3$color
names(colors3) <- rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq)]))


## We create the plot
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_plot <- ggplot(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#370335B2",
                               "Tomato-Ripening" = "#bea3ce",#"#bebada",
                               "Pepper-Ripening" = "#fccde5",
                               "Ripening"="#bc80bd",
                               "Nicotiana-independent"="#808080",#"#D2AF81B2",
                               "OR-Dependent"="#F05C3BB2", 
                               "Ripening-independent"="#FED439B2",
                               "Fruit-independent"="#006c00",
                               "Spread"="#b3de69", 
                               "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                               "Exclusive"="#197EC0B2"),
                     breaks = c("Constant",
                                "Tomato-Ripening",
                                "Pepper-Ripening",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent",
                                "Fruit-independent",
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Tomato-Pepper-Nicotiana-Arabidopsis Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme(axis.text.y = element_text(color = colors3, face = "bold"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_plot ## Saved pdf 100x28 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Second approach

for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up)){
  if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]==1){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Exclusive"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & 
            any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(1, 4)) & 
            #any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(7)) & 
            #any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(9)) & 
            all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %ni% c(2,3,5,6)) & 
            any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(11:18)) &
            all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %ni% c(19))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Early"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#9d8545"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$NumberSets[i]>1 & 
            all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %ni% c(1, 4)) & 
            any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(2, 3, 5, 6)) & 
            #any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(8)) & 
            #any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(10)) & 
            all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %ni% c(11:18)) & 
            any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$TotalSets[i])) %in% c(19))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Late"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#d95f02"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq[i] <- "Whatever"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq, levels = c("Early", "Late","Exclusive", "Whatever"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description, levels = (rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors3 <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up %>% dplyr::select(description, color)
colors3 <- unique(colors3)
colors3 <- colors3[order(colors3$description),]
colors3 <- colors3$color
names(colors3) <- rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$Group, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq)]))


## We create the plot
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_plot <- ggplot(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up %>%
                                                                                   filter(freq %in% c("Early", "Late")), aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Early"="#9d8545",
                               "Late"="#d95f02"),
                     breaks = c("Early",
                                "Late"))+
  ggtitle("Biological Processes Tomato-Pepper-Nicotiana-Arabidopsis Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme(axis.text.y = element_text(color = colors3, face = "bold"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_plot ## Saved pdf 100x28 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


########### Down

filelist_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down <- lapply(list("topGO_weightFS_Pk_Down_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_LR_Down_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_RR_Down_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_M82_B_Down_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_M82_or_Down_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_M82_Red_Down_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_Breaker_Down_BP_REVIGO",
                                                                          "topGO_weightFS_Red_Down_BP_REVIGO",
                                                                          "topGO_weightFS_MG_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_B_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                                          "topGO_weightFS_22h_Down_BP_REVIGO",
                                                                          "topGO_weightFS_25h_Down_BP_REVIGO",
                                                                          "topGO_weightFS_28h_Down_BP_REVIGO",
                                                                          "topGO_weightFS_34h_Down_BP_REVIGO",
                                                                          "topGO_weightFS_37h_Down_BP_REVIGO",
                                                                          "topGO_weightFS_40h_Down_BP_REVIGO",
                                                                          "topGO_weightFS_46h_Down_BP_REVIGO",
                                                                          "topGO_weightFS_56h_Down_BP_REVIGO",
                                                                          "topGO_weightFS_96h_Down_BP_REVIGO",
                                                                          "topGO_weightFS_Woo_Down_BP_REVIGO"), get)
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down <- do.call(rbind, filelist_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down)

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Stage <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Stage,
                                                                                  levels = c("Pk", "LR", "RR", "B", "or", "Red", "Breaker", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h", "Woo"),
                                                                                  labels = c("Breaker", "Orange", "Red", "Breaker", "Orange", "Red", "Breaker", "Red", "ORHis", "22h", "25h", "28h","34h","37h","40h","46h","56h","96h", "Senescence"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Background <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Background,
                                                                                       levels = c("M82", "Pepper","MG", "B", "Niben", "Ath"),
                                                                                       labels = c("M82", "Pepper","MG", "Breaker", "Niben", "Ath"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group <- paste(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Stage, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Background, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Set, sep = "_")
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group,
                                                                                  levels = c("Breaker_M82_Shinozaki", "Orange_M82_Shinozaki", "Red_M82_Shinozaki",
                                                                                             "Breaker_M82_Yazdani", "Orange_M82_Yazdani", "Red_M82_Yazdani",
                                                                                             "Breaker_Pepper_Capsicum", "Red_Pepper_Capsicum","ORHis_MG_Yazdani", "ORHis_Breaker_Yazdani",
                                                                                             "22h_Niben_Nicotiana", "25h_Niben_Nicotiana", "28h_Niben_Nicotiana", "34h_Niben_Nicotiana", 
                                                                                             "37h_Niben_Nicotiana", "40h_Niben_Nicotiana", "46h_Niben_Nicotiana", "56h_Niben_Nicotiana", 
                                                                                             "96h_Niben_Nicotiana", "Senescence_Ath_Arabidopsis"),
                                                                                  labels = c("Breaker_Shinozaki", "Orange_Shinozaki", "Red_Shinozaki",
                                                                                             "Breaker_Yazdani", "Orange_Yazdani", "Red_Yazdani",
                                                                                             "Breaker_Pepper", "Red_Pepper","ORHis_MG", "ORHis_Breaker",
                                                                                             "22h_Niben", "25h_Niben", "28h_Niben", "34h_Niben", 
                                                                                             "37h_Niben", "40h_Niben", "46h_Niben", "56h_Niben", 
                                                                                             "96h_Niben", "Senescence_Ath"))

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group),]

for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down)){
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i] <- list(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group[Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description==Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description[i]])
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i] <- length(c(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group[Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description==Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description[i]]))
}


# First approach
for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down)){
  if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]==1){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Exclusive"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>13){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Constant"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#370335B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))<7)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Tomato-Ripening"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#bea3ce"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))>6) & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))<9)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Pepper-Ripening"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#fccde5"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))<9) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(1:6)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(7,8))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Ripening"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#bc80bd"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))>8) & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))<11)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "OR-Dependent"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#F05C3BB2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))>10) & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))<20)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Nicotiana"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#D2AF81B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))>8) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(9,10)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))>10)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Ripening-independent"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#FED439B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(1:10)) & all(!as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(11:19))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#808080"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))>10) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(11:19)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i]))>19)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Fruit-independent"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#006c00"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Spread"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq, levels = c("Constant","Tomato-Ripening", "Pepper-Ripening", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Fruit-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description, levels = (rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors4 <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down %>% dplyr::select(description, color)
colors4 <- unique(colors4)
colors4 <- colors4[order(colors4$description),]
colors4 <- colors4$color
names(colors4) <- rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq)]))


## We create the plot
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_plot <- ggplot(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#370335B2",
                               "Tomato-Ripening" = "#bea3ce",#"#bebada",
                               "Pepper-Ripening" = "#fccde5",
                               "Ripening"="#bc80bd",
                               "Nicotiana-independent"="#808080",#"#D2AF81B2",
                               "OR-Dependent"="#F05C3BB2", 
                               "Ripening-independent"="#FED439B2",
                               "Fruit-independent"="#006c00",
                               "Spread"="#b3de69", 
                               "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                               "Exclusive"="#197EC0B2"),
                     breaks = c("Constant",
                                "Tomato-Ripening",
                                "Pepper-Ripening",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent",
                                "Fruit-independent",
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Tomato-Pepper-Nicotiana-Arabidopsis Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme(axis.text.y = element_text(color = colors4, face = "bold"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_plot ## Saved pdf 100x28 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


# Second approach

for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down)){
  if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]==1){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Exclusive"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & 
            any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(1, 4)) & 
            #any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(7)) & 
            #any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(9)) & 
            all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %ni% c(2,3,5,6)) & 
            any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(11:18)) &
            all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %ni% c(19))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Early"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#9d8545"
  } else if(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$NumberSets[i]>1 & 
            all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %ni% c(1, 4)) & 
            any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(2, 3, 5, 6)) & 
            #any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(8)) & 
            #any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(10)) & 
            all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %ni% c(11:18)) & 
            any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$TotalSets[i])) %in% c(19))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Late"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#d95f02"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq[i] <- "Whatever"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq, levels = c("Early", "Late","Exclusive", "Whatever"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description, levels = (rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq)]))))

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down %>%
  filter(freq %in% c("Early", "Late"))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors3 <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down %>% dplyr::select(description, color)
colors3 <- unique(colors3)
colors3 <- colors3[order(colors3$description),]
colors3 <- colors3$color
names(colors3) <- rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$Group, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq)]))


## We create the plot
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_plot <- ggplot(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down %>%
                                                                                   filter(freq %in% c("Early", "Late")), aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Early"="#9d8545",
                               "Late"="#d95f02"),
                     breaks = c("Early",
                                "Late"))+
  ggtitle("Biological Processes Tomato-Pepper-Nicotiana-Arabidopsis Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme(axis.text.y = element_text(color = colors3, face = "bold"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_plot ## Saved pdf 100x28 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)




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

######################################################
## SUMMARIZING RESULTS BY GROUP (% AND COUNTS)

## Up
percentData_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up %>%
  group_by(Group) %>% count(freq) %>%
  mutate(ratio=scales::percent(n/sum(n)))
percentData_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_percent <- ggplot(
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up, aes(x= Group))+
  geom_bar(aes(fill=freq), position = position_fill(reverse = TRUE))+
  geom_text(data = percentData_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up, aes(y=n, label=ratio),
            position = position_fill(vjust = 0.5), colour = "white", fontface = "bold")+
  scale_fill_manual(values= c("Constant"="#370335B2",
                              "Tomato-Ripening" = "#bea3ce",#"#bebada",
                              "Pepper-Ripening" = "#fccde5",
                              "Ripening"="#bc80bd",
                              "Nicotiana-independent"="#808080",#"#D2AF81B2",
                              "OR-Dependent"="#F05C3BB2", 
                              "Ripening-independent"="#FED439B2",
                              "Fruit-independent"="#006c00",
                              "Spread"="#b3de69", 
                              "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                              "Exclusive"="#197EC0B2"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     labels = c("0", "25", "50", "75", "100"))+
  ylab("% of GO terms subgroups")+
  theme_HPLC+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_percent 

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_abs <- ggplot(
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up, aes(x=Group))+
  geom_bar(aes(fill=Background))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_discrete(labels = c("Breaker_Shino", "Orange_Shino", "Red_Shino","Breaker_Yaz", "Orange_Yaz", "Red_Yaz", "Breaker_Pep", "Red_Pep", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h", "Ath"))+
  ylab("Number of GO terms")+
  scale_fill_manual(values = c("M82" = "#bea3ce",
                               "Pepper" ="#fccde5",
                               "MG"="#e45132",
                               "Breaker"="#bd2e15",
                               "Niben"="#D2AF81B2",
                               "Ath"= "#006c00"))+
  theme_HPLC+
  theme(legend.position = "top",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_abs

figureA <- ggarrange(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_abs,
                     Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_percent,
                     ncol = 1,
                     nrow = 2)
figureA ## Saved pdf 15x10 "Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Up_Percent"


## Down
percentData_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down %>%
  group_by(Group) %>% count(freq) %>%
  mutate(ratio=scales::percent(n/sum(n)))
percentData_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_percent <- ggplot(
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down, aes(x= Group))+
  geom_bar(aes(fill=freq), position = position_fill(reverse = TRUE))+
  geom_text(data = percentData_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down, aes(y=n, label=ratio),
            position = position_fill(vjust = 0.5), colour = "white", fontface = "bold")+
  scale_fill_manual(values= c("Constant"="#370335B2",
                              "Tomato-Ripening" = "#bea3ce",#"#bebada",
                              "Pepper-Ripening" = "#fccde5",
                              "Ripening"="#bc80bd",
                              "Nicotiana-independent"="#808080",#"#D2AF81B2",
                              "OR-Dependent"="#F05C3BB2", 
                              "Ripening-independent"="#FED439B2",
                              "Fruit-independent"="#006c00",
                              "Spread"="#b3de69", 
                              "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                              "Exclusive"="#197EC0B2"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     labels = c("0", "25", "50", "75", "100"))+
  ylab("% of GO terms subGroups")+
  theme_HPLC+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_percent 

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_abs <- ggplot(
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down, aes(x=Group))+
  geom_bar(aes(fill=Background))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_discrete(labels = c("Breaker_Shino", "Orange_Shino", "Red_Shino","Breaker_Yaz", "Orange_Yaz", "Red_Yaz", "Breaker_Pep", "Red_Pep", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h", "Ath"))+
  ylab("Number of GO terms")+
  scale_fill_manual(values = c("M82" = "#bea3ce",
                               "Pepper" ="#fccde5",
                               "MG"="#e45132",
                               "Breaker"="#bd2e15",
                               "Niben"="#D2AF81B2",
                               "Ath"= "#006c00"))+
  theme_HPLC+
  theme(legend.position = "top",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_abs

figureB <- ggarrange(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_abs,
                     Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_percent,
                     ncol = 1,
                     nrow = 2)
figureB ## Saved pdf 15x10 "Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_Down_Percent"


############################################################################
## Graph only Spread

## Up

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up[Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq=="Spread",]

'%ni%' <- Negate('%in%')
# Try 1
for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i]))<16) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i])) %in% c(1:8)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i])) %in% c(11:15))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2[i] <- "Early"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i]))%in% c(1:8)) & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i])) %ni% c(11:15)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i])) %in% c(16:19))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2[i] <- "Late"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$color2[i] <- "#7570b3"}
}

#Try 2

for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i]))!=19)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2[i] <- "Early"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$color2[i] <- "#1b9e77"
  } else if(all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i]))%ni% c(11:18)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i]))== 19)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2[i] <- "Late"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$color2[i] <- "#7570b3"}
}


Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2 <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$description <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$description, levels = (rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$Group)]))))

colors5 <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread %>% dplyr::select(description, color2)
colors5 <- unique(colors5)
colors5 <- colors5[order(colors5$description),]
colors5 <- colors5$color2
names(colors5) <- rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$Group)]))


Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_plot <- ggplot(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 8.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 8.5, xmax = 9.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 9.5, xmax = 10.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 10.5, xmax = 19.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  annotate("rect",xmin = 19.5, xmax = 20.5,
           ymin = -Inf, ymax = Inf,
           fill = "#006c00", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,20,1),
                     labels = c("Breaker_Shino", "Orange_Shino", "Red_Shino","Breaker_Yaz", "Orange_Yaz", "Red_Yaz", "Breaker_Pep", "Red_Pep", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h", "Ath"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Shinozaki-Yazdani-Pepper-Nicotiana-Ath Spread Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors5, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_plot ## Saved pdf 40x20 "Summary_topGO_Shinzaki&Yazdani&Pepper&Nb&Ath_Up_0.7_Spread"

write.table((Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread %>% dplyr::select(description, Group, freq2)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Shinzaki&Yazdani&Pepper&Nb&Ath_Up_0.7_Spread.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(unique((Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread %>% dplyr::select(description, freq2) %>% dplyr::arrange(freq2, description))),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Shinzaki&Yazdani&Pepper&Nb&Ath_Up_0.7_Spread_2.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

##Alternative

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up[Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up$freq=="Spread",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$TotalSets[i]))!=19)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2[i] <- "Early"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$color2[i] <- "#1b9e77"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$freq2[i] <- "Late"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread$color2[i] <- "#d95f02"
  }
}

Selection_Up <- c("protein refolding",
                  "oxidation-reduction process",
                  "chorismate biosynthetic process",
                  "chaperone mediated protein folding requiring cofactor"
)

Early_Up <- c("regulation of sulfur utilization",
              "cellular response to sulfur starvation",
              "positive regulation of response to water deprivation",
              "dimethylallyl diphosphate biosynthetic process",
              "piecemeal microautophagy of nucleus")
Late_Up <- c("positive regulation of necrotic cell death",
             "negative regulation of mitochondrial membrane potential",
             "protein stabilization",
             "'de novo' IMP biosynthetic process",
             "chloroplast rRNA processing",
             "'de novo' CTP biosynthetic process",
             "chloroplast mRNA processing",
             "protein maturation by iron-sulfur cluster transfer",
             "response to heat",
             "positive regulation of superoxide dismutase activity")

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2 <-
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread %>%
  filter(description %in% c(Early_Up, Late_Up))

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$freq2 <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$freq2, levels = c("Early", "Late"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$description <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$description, levels = (rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$freq2, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$Group)]))))

colors5 <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2 %>% dplyr::select(description, color2)
colors5 <- unique(colors5)
colors5 <- colors5[order(colors5$description),]
colors5 <- colors5$color2
names(colors5) <- rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$freq2, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2$Group)]))



Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_plot <- ggplot(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_2, 
                                                                                        aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 8.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 8.5, xmax = 9.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 9.5, xmax = 10.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 10.5, xmax = 19.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  annotate("rect",xmin = 19.5, xmax = 20.5,
           ymin = -Inf, ymax = Inf,
           fill = "#006c00", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02"),
                     breaks = c("Early",
                                "Late"))+
  scale_x_continuous(breaks = seq(1,20,1),
                     labels = c("Breaker_Shino", "Orange_Shino", "Red_Shino","Breaker_Yaz", "Orange_Yaz", "Red_Yaz", "Breaker_Pep", "Red_Pep", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h", "Ath"),
                     expand = expansion(mult = 0, add = 0))+
  #ggtitle("Biological Processes Shinozaki-Yazdani-Pepper-Nicotiana-Ath Spread Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors5, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank(),
        axis.text.y.left = element_blank(),
        legend.position = "bottom")
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Up_Spread_plot ## Saved pdf 40x20 "Summary_topGO_Shinzaki&Yazdani&Pepper&Nb&Ath_Up_0.7_Spread"


## Down

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down[Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq=="Spread",]

'%ni%' <- Negate('%in%')

# Try 1
for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i]))<16) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i])) %in% c(1:8)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i])) %in% c(11:15))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2[i] <- "Early"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i]))%in% c(1:8)) & all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i])) %ni% c(11:15)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i])) %in% c(16:19))){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2[i] <- "Late"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$color2[i] <- "#7570b3"}
}

#Try 2

for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i]))!=19)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2[i] <- "Early"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$color2[i] <- "#1b9e77"
  } else if(all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i]))%ni% c(11:18)) & any(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i]))== 19)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2[i] <- "Late"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$color2[i] <- "#7570b3"}
}

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2 <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$description <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$description, levels = (rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$Group)]))))

colors6 <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread %>% dplyr::select(description, color2)
colors6 <- unique(colors6)
colors6 <- colors6[order(colors6$description),]
colors6 <- colors6$color2
names(colors6) <- rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$Group)]))


Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_plot <- ggplot(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 8.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 8.5, xmax = 9.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 9.5, xmax = 10.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 10.5, xmax = 19.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  annotate("rect",xmin = 19.5, xmax = 20.5,
           ymin = -Inf, ymax = Inf,
           fill = "#006c00", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,20,1),
                     labels = c("Breaker_Shino", "Orange_Shino", "Red_Shino","Breaker_Yaz", "Orange_Yaz", "Red_Yaz", "Breaker_Pep", "Red_Pep", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h", "Ath"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Shinozaki-Yazdani-Pepper-Nicotiana-Ath Spread Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors6, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_plot ## Saved pdf 30x20 "Summary_topGO_Shinzaki&Yazdani&Pepper&Nb&Ath_Up_0.7_Spread"

write.table((Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread %>% dplyr::select(description, Group, freq2)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Shinzaki&Yazdani&Pepper&Nb&Ath_Down_0.7_Spread.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(unique((Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread %>% dplyr::select(description, freq2) %>% dplyr::arrange(freq2, description))),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Shinzaki&Yazdani&Pepper&Nb&Ath_Down_0.7_Spread_2.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)



##Alternative

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down[Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down$freq=="Spread",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$TotalSets[i]))<19)){
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2[i] <- "Early"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$color2[i] <- "#1b9e77"
  } else {
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$freq2[i] <- "Late"
    Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread$color2[i] <- "#d95f02"
  }
}

Selection_Down <- c("anthocyanin-containing compound biosynthetic process",
                    "transmembrane transport",
                    "regulation of transcription, DNA-templated",
                    "pre-replicative complex assembly involved in nuclear cell cycle DNA replication"
)

Early_Down <- c("trehalose biosynthetic process",
                "L-serine biosynthetic process",
                "ion transmembrane transport",
                "vegetative phase change",
                "pre-replicative complex assembly involved in nuclear cell cycle DNA replication")
Late_Down <- c("male meiosis II",
               "cytokinin-activated signaling pathway",
               "sucrose transport",
               "cellulose biosynthetic process")

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2 <-
  Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread %>%
  filter(description %in% c(Early_Down, Late_Down))

Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$freq2 <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$freq2, levels = c("Early", "Late"))
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$description <- factor(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$description, levels = (rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$freq2, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$Group)]))))

colors5 <- Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2 %>% dplyr::select(description, color2)
colors5 <- unique(colors5)
colors5 <- colors5[order(colors5$description),]
colors5 <- colors5$color2
names(colors5) <- rev(unique(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$description[order(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$freq2, Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2$Group)]))



Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_plot <- ggplot(Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_2, 
                                                                                        aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 8.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 8.5, xmax = 9.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 9.5, xmax = 10.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 10.5, xmax = 19.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  annotate("rect",xmin = 19.5, xmax = 20.5,
           ymin = -Inf, ymax = Inf,
           fill = "#006c00", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02"),
                     breaks = c("Early",
                                "Late"))+
  scale_x_continuous(breaks = seq(1,20,1),
                     labels = c("Breaker_Shino", "Orange_Shino", "Red_Shino","Breaker_Yaz", "Orange_Yaz", "Red_Yaz", "Breaker_Pep", "Red_Pep", "MG", "Breaker", "22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h", "Ath"),
                     expand = expansion(mult = 0, add = 0))+
  #ggtitle("Biological Processes Shinozaki-Yazdani-Pepper-Nicotiana-Ath Spread Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors5, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank(),
        axis.text.y.left = element_blank(),
        legend.position = "none")
Summary_GOterms_Shinozaki_vs_Yazdani_vs_Pepper_vs_Niben_vs_Ath_Down_Spread_plot ## Saved pdf 40x20 "Summary_topGO_Shinzaki&Yazdani&Pepper&Nb&Ath_Down_0.7_Spread"








##################################################################################################
########################################  PCA  ###########################################
##################################################################################################
## YAZDANI VS NICOTIANA

########### Up

filelist_Yazdani_vs_Niben_PCA_Up <- lapply(list("topGO_weightFS_M82_B_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_M82_or_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_M82_Red_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_MG_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_B_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                            "topGO_weightFS_AA_Up_BP_REVIGO",
                                            "topGO_weightFS_BB_Up_BP_REVIGO",
                                            "topGO_weightFS_CC_Up_BP_REVIGO",
                                            "topGO_weightFS_46h_Up_BP_REVIGO",
                                            "topGO_weightFS_56h_Up_BP_REVIGO",
                                            "topGO_weightFS_96h_Up_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Niben_PCA_Up <- do.call(rbind, filelist_Yazdani_vs_Niben_PCA_Up)

Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Stage <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Stage,
                                                    levels = c("B", "or", "Red", "ORHis", "AA", "BB", "CC","46h","56h","96h"),
                                                    labels = c("Breaker", "Orange", "Red", "ORHis", "AA", "BB", "CC","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Background <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Background,
                                                         levels = c("M82", "MG", "B", "Niben"),
                                                         labels = c("M82", "MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Group <- paste(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Stage, Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Group <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Group,
                                                    levels = c("Breaker_M82", "Orange_M82", "Red_M82", "ORHis_MG", "ORHis_Breaker","AA_Niben", "BB_Niben", "CC_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Niben_PCA_Up <- Summary_GOterms_Yazdani_vs_Niben_PCA_Up[order(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_PCA_Up)){
  Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Group[Summary_GOterms_Yazdani_vs_Niben_PCA_Up$description==Summary_GOterms_Yazdani_vs_Niben_PCA_Up$description[i]])
  Summary_GOterms_Yazdani_vs_Niben_PCA_Up$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Group[Summary_GOterms_Yazdani_vs_Niben_PCA_Up$description==Summary_GOterms_Yazdani_vs_Niben_PCA_Up$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_PCA_Up)){
  if(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$color[i] <- "#8DA0CB"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$NumberSets[i]>8){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$color[i] <- "#66C2A5"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$color[i] <- "#E78AC3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$color[i] <- "#FC8D62"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i]))>5)){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$color[i] <- "#B3B3B3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i]))>3) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i])) %in% c(6:11))){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$color[i] <- "#FFD92F"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$TotalSets[i])) %in% c(1:3))){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$color[i] <- "#E5C494"
  } else {
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Up$color[i] <- "#A6D854"}
}

## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq, levels = c("Constant", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Niben_PCA_Up$description <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$description[order(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Group, Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors9 <- Summary_GOterms_Yazdani_vs_Niben_PCA_Up %>% dplyr::select(description, color)
colors9 <- unique(colors9)
colors9 <- colors9[order(colors9$description),]
colors9 <- colors9$color
names(colors9) <- rev(unique(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$description[order(Summary_GOterms_Yazdani_vs_Niben_PCA_Up$Group, Summary_GOterms_Yazdani_vs_Niben_PCA_Up$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Niben_PCA_Up_plot <- ggplot(Summary_GOterms_Yazdani_vs_Niben_PCA_Up, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#66C2A5",
                               "Ripening"="#E78AC3",
                               "Nicotiana-independent"="#E5C494",
                               "OR-Dependent"="#FC8D62", 
                               "Ripening-independent"="#FFD92F", 
                               "Spread"="#A6D854", 
                               "Nicotiana"="#B3B3B3",
                               "Exclusive"="#8DA0CB"),
                     breaks = c("Constant",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Nicotiana Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme(axis.text.y = element_text(color = colors9, face = "bold"))
Summary_GOterms_Yazdani_vs_Niben_PCA_Up_plot ## Saved pdf 20x70 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Yazdani_vs_Niben_PCA_Up %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Yazdani_vs_Niben_PCA_Up_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


########### Down

filelist_Yazdani_vs_Niben_PCA_Down <- lapply(list("topGO_weightFS_M82_B_Down_BP_4.0_AHRD_REVIGO",
                                              "topGO_weightFS_M82_or_Down_BP_4.0_AHRD_REVIGO",
                                              "topGO_weightFS_M82_Red_Down_BP_4.0_AHRD_REVIGO",
                                              "topGO_weightFS_MG_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                              "topGO_weightFS_B_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                              "topGO_weightFS_AA_Down_BP_REVIGO",
                                              "topGO_weightFS_BB_Down_BP_REVIGO",
                                              "topGO_weightFS_CC_Down_BP_REVIGO",
                                              "topGO_weightFS_46h_Down_BP_REVIGO",
                                              "topGO_weightFS_56h_Down_BP_REVIGO",
                                              "topGO_weightFS_96h_Down_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Niben_PCA_Down <- do.call(rbind, filelist_Yazdani_vs_Niben_PCA_Down)

Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Stage <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Stage,
                                                      levels = c("B", "or", "Red", "ORHis", "AA", "BB", "CC","46h","56h","96h"),
                                                      labels = c("Breaker", "Orange", "Red", "ORHis", "AA", "BB", "CC","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Background <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Background,
                                                           levels = c("M82", "MG", "B", "Niben"),
                                                           labels = c("M82", "MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Group <- paste(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Stage, Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Group <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Group,
                                                      levels = c("Breaker_M82", "Orange_M82", "Red_M82", "ORHis_MG", "ORHis_Breaker","AA_Niben", "BB_Niben", "CC_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Niben_PCA_Down <- Summary_GOterms_Yazdani_vs_Niben_PCA_Down[order(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_PCA_Down)){
  Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Group[Summary_GOterms_Yazdani_vs_Niben_PCA_Down$description==Summary_GOterms_Yazdani_vs_Niben_PCA_Down$description[i]])
  Summary_GOterms_Yazdani_vs_Niben_PCA_Down$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Group[Summary_GOterms_Yazdani_vs_Niben_PCA_Down$description==Summary_GOterms_Yazdani_vs_Niben_PCA_Down$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_PCA_Down)){
  if(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$color[i] <- "#8DA0CB"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$NumberSets[i]>8){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$color[i] <- "#66C2A5"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$color[i] <- "#E78AC3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$color[i] <- "#FC8D62"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i]))>5)){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$color[i] <- "#B3B3B3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i]))>3) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i])) %in% c(6:11))){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$color[i] <- "#FFD92F"
  } else if(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$TotalSets[i])) %in% c(1:3))){
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$color[i] <- "#E5C494"
  } else {
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Niben_PCA_Down$color[i] <- "#A6D854"}
}

## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq, levels = c("Constant", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Niben_PCA_Down$description <- factor(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$description[order(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Group, Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors10 <- Summary_GOterms_Yazdani_vs_Niben_PCA_Down %>% dplyr::select(description, color)
colors10 <- unique(colors10)
colors10 <- colors10[order(colors10$description),]
colors10 <- colors10$color
names(colors10) <- rev(unique(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$description[order(Summary_GOterms_Yazdani_vs_Niben_PCA_Down$Group, Summary_GOterms_Yazdani_vs_Niben_PCA_Down$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Niben_PCA_Down_plot <- ggplot(Summary_GOterms_Yazdani_vs_Niben_PCA_Down, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#66C2A5",
                               "Ripening"="#E78AC3",
                               "Nicotiana-independent"="#E5C494",
                               "OR-Dependent"="#FC8D62", 
                               "Ripening-independent"="#FFD92F", 
                               "Spread"="#A6D854", 
                               "Nicotiana"="#B3B3B3",
                               "Exclusive"="#8DA0CB"),
                     breaks = c("Constant",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Nicotiana Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme(axis.text.y = element_text(color = colors10, face = "bold"))
Summary_GOterms_Yazdani_vs_Niben_PCA_Down_plot ## Saved pdf 20x70 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Yazdani_vs_Niben_PCA_Down %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Yazdani_vs_Niben_PCA_Down_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)



##################################################################################################
## YAZDANI VS PEPPER VS NICOTIANA

########### Up

filelist_Yazdani_vs_Pepper_vs_Niben_PCA_Up <- lapply(list("topGO_weightFS_M82_B_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_M82_or_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_M82_Red_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_Breaker_Up_BP_REVIGO",
                                                      "topGO_weightFS_Red_Up_BP_REVIGO",
                                                      "topGO_weightFS_MG_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_B_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                      "topGO_weightFS_AA_Up_BP_REVIGO",
                                                      "topGO_weightFS_BB_Up_BP_REVIGO",
                                                      "topGO_weightFS_CC_Up_BP_REVIGO",
                                                      "topGO_weightFS_46h_Up_BP_REVIGO",
                                                      "topGO_weightFS_56h_Up_BP_REVIGO",
                                                      "topGO_weightFS_96h_Up_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up <- do.call(rbind, filelist_Yazdani_vs_Pepper_vs_Niben_PCA_Up)

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Stage <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Stage,
                                                              levels = c("B", "or", "Red", "Breaker", "Red", "ORHis", "AA", "BB", "CC","46h","56h","96h"),
                                                              labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "ORHis", "AA", "BB", "CC","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Background <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Background,
                                                                   levels = c("M82", "Pepper","MG", "B", "Niben"),
                                                                   labels = c("M82", "Pepper","MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Group <- paste(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Stage, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Group <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Group,
                                                              levels = c("Breaker_M82", "Orange_M82", "Red_M82", "Breaker_Pepper", "Red_Pepper","ORHis_MG", "ORHis_Breaker","AA_Niben", "BB_Niben", "CC_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up)){
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$description[i]])
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up)){
  if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i]>10){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#370335B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "Tomato-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#bea3ce"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "Pepper-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#fccde5"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i])) %in% c(1:3)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i])) %in% c(4,5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#bc80bd"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i]))>5) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i]))<8)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#F05C3BB2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i]))>7)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#D2AF81B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i]))>5) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i])) %in% c(8:13))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#FED439B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i]))<8) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$TotalSets[i])) %in% c(1:5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#808080"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq, levels = c("Constant","Tomato-Ripening", "Pepper-Ripening", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors11 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up %>% dplyr::select(description, color)
colors11 <- unique(colors11)
colors11 <- colors11[order(colors11$description),]
colors11 <- colors11$color
names(colors11) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#370335B2",
                               "Tomato-Ripening" = "#bea3ce",#"#bebada",
                               "Pepper-Ripening" = "#fccde5",
                               "Ripening"="#bc80bd",
                               "Nicotiana-independent"="#808080",#"#D2AF81B2",
                               "OR-Dependent"="#F05C3BB2", 
                               "Ripening-independent"="#FED439B2", 
                               "Spread"="#b3de69", 
                               "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                               "Exclusive"="#197EC0B2"),
                     breaks = c("Constant",
                                "Tomato-Ripening",
                                "Pepper-Ripening",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme(axis.text.y = element_text(color = colors11, face = "bold"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_plot ## Saved pdf 20x73 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

########### Down

filelist_Yazdani_vs_Pepper_vs_Niben_PCA_Down <- lapply(list("topGO_weightFS_M82_B_Down_BP_4.0_AHRD_REVIGO",
                                                        "topGO_weightFS_M82_or_Down_BP_4.0_AHRD_REVIGO",
                                                        "topGO_weightFS_M82_Red_Down_BP_4.0_AHRD_REVIGO",
                                                        "topGO_weightFS_Breaker_Down_BP_REVIGO",
                                                        "topGO_weightFS_Red_Down_BP_REVIGO",
                                                        "topGO_weightFS_MG_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                        "topGO_weightFS_B_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                        "topGO_weightFS_AA_Down_BP_REVIGO",
                                                        "topGO_weightFS_BB_Down_BP_REVIGO",
                                                        "topGO_weightFS_CC_Down_BP_REVIGO",
                                                        "topGO_weightFS_46h_Down_BP_REVIGO",
                                                        "topGO_weightFS_56h_Down_BP_REVIGO",
                                                        "topGO_weightFS_96h_Down_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down <- do.call(rbind, filelist_Yazdani_vs_Pepper_vs_Niben_PCA_Down)

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Stage <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Stage,
                                                                levels = c("B", "or", "Red", "Breaker", "Red", "ORHis", "AA", "BB", "CC","46h","56h","96h"),
                                                                labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "ORHis", "AA", "BB", "CC","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Background <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Background,
                                                                     levels = c("M82", "Pepper","MG", "B", "Niben"),
                                                                     labels = c("M82", "Pepper","MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Group <- paste(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Stage, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Group <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Group,
                                                                levels = c("Breaker_M82", "Orange_M82", "Red_M82", "Breaker_Pepper", "Red_Pepper","ORHis_MG", "ORHis_Breaker","AA_Niben", "BB_Niben", "CC_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down)){
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$description[i]])
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down)){
  if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i]>10){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#370335B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "Tomato-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#bea3ce"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "Pepper-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#fccde5"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i])) %in% c(1:3)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i])) %in% c(4,5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#bc80bd"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i]))>5) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i]))<8)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#F05C3BB2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i]))>7)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#D2AF81B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i]))>5) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i])) %in% c(8:13))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#FED439B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i]))<8) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$TotalSets[i])) %in% c(1:5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#808080"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq, levels = c("Constant","Tomato-Ripening", "Pepper-Ripening", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors12 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down %>% dplyr::select(description, color)
colors12 <- unique(colors12)
colors12 <- colors12[order(colors12$description),]
colors12 <- colors12$color
names(colors12) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#370335B2",
                               "Tomato-Ripening" = "#bea3ce",#"#bebada",
                               "Pepper-Ripening" = "#fccde5",
                               "Ripening"="#bc80bd",
                               "Nicotiana-independent"="#808080",#"#D2AF81B2",
                               "OR-Dependent"="#F05C3BB2", 
                               "Ripening-independent"="#FED439B2", 
                               "Spread"="#b3de69", 
                               "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                               "Exclusive"="#197EC0B2"),
                     breaks = c("Constant",
                                "Tomato-Ripening",
                                "Pepper-Ripening",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme(axis.text.y = element_text(color = colors12, face = "bold"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_plot ## Saved pdf 20x73 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

######################################################
## SUMMARIZING RESULTS BY GROUP (% AND COUNTS)

## Up
percentData_Yazdani_vs_Pepper_vs_Niben_PCA_Up <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up %>%
  group_by(Group) %>% count(freq) %>%
  mutate(ratio=scales::percent(n/sum(n)))
percentData_Yazdani_vs_Pepper_vs_Niben_PCA_Up

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_percent <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up, aes(x= Group))+
  geom_bar(aes(fill=freq), position = position_fill(reverse = TRUE))+
  geom_text(data = percentData_Yazdani_vs_Pepper_vs_Niben_PCA_Up, aes(y=n, label=ratio),
            position = position_fill(vjust = 0.5), colour = "white", fontface = "bold")+
  scale_fill_manual(values= c("Constant"="#370335B2",
                              "Tomato-Ripening" = "#bea3ce",#"#bebada",
                              "Pepper-Ripening" = "#fccde5",
                              "Ripening"="#bc80bd",
                              "Nicotiana-independent"="#808080",#"#D2AF81B2",
                              "OR-Dependent"="#F05C3BB2", 
                              "Ripening-independent"="#FED439B2", 
                              "Spread"="#b3de69", 
                              "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                              "Exclusive"="#197EC0B2"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     labels = c("0", "25", "50", "75", "100"))+
  ylab("% of GO terms subgroups")+
  theme_HPLC+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_percent 

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_abs <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up, aes(x=Group))+
  geom_bar(aes(fill=Background))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_discrete(labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "AA", "BB", "CC", "46h", "56h", "96h"))+
  ylab("Number of GO terms")+
  scale_fill_manual(values = c("M82" = "#bea3ce",
                               "Pepper" ="#fccde5",
                               "MG"="#e45132",
                               "Breaker"="#bd2e15",
                               "Niben"="#D2AF81B2"))+
  theme_HPLC+
  theme(legend.position = "top",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_abs

figureC <- ggarrange(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_abs,
                     Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_percent,
                     ncol = 1,
                     nrow = 2)
figureC ## Saved pdf 15x10 "Summary_topGO_Yazdani&Pepper&NbPCA_Up_0.7_Percent"


## Down
percentData_Yazdani_vs_Pepper_vs_Niben_PCA_Down <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down %>%
  group_by(Group) %>% count(freq) %>%
  mutate(ratio=scales::percent(n/sum(n)))
percentData_Yazdani_vs_Pepper_vs_Niben_PCA_Down

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_percent <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down, aes(x= Group))+
  geom_bar(aes(fill=freq), position = position_fill(reverse = TRUE))+
  geom_text(data = percentData_Yazdani_vs_Pepper_vs_Niben_PCA_Down, aes(y=n, label=ratio),
            position = position_fill(vjust = 0.5), colour = "white", fontface = "bold")+
  scale_fill_manual(values= c("Constant"="#370335B2",
                              "Tomato-Ripening" = "#bea3ce",#"#bebada",
                              "Pepper-Ripening" = "#fccde5",
                              "Ripening"="#bc80bd",
                              "Nicotiana-independent"="#808080",#"#D2AF81B2",
                              "OR-Dependent"="#F05C3BB2", 
                              "Ripening-independent"="#FED439B2", 
                              "Spread"="#b3de69", 
                              "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                              "Exclusive"="#197EC0B2"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     labels = c("0", "25", "50", "75", "100"))+
  ylab("% of GO terms subGroups")+
  theme_HPLC+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_percent 

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_abs <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down, aes(x=Group))+
  geom_bar(aes(fill=Background))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_discrete(labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "AA", "BB", "CC", "46h", "56h", "96h"))+
  scale_fill_manual(values = c("M82" = "#bea3ce",
                               "Pepper" ="#fccde5",
                               "MG"="#e45132",
                               "Breaker"="#bd2e15",
                               "Niben"="#D2AF81B2"))+
  theme_HPLC+
  theme(legend.position = "top",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_abs

figureD <- ggarrange(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_abs,
                     Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_percent,
                     ncol = 1,
                     nrow = 2)
figureD ## Saved pdf 15x10 "Summary_topGO_Yazdani&Pepper&NbPCA_Down_0.7_Percent"



############################################################################
## Graph only Spread

## Up

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq=="Spread",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$TotalSets[i]))<11) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$TotalSets[i])) %in% c(1:5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$TotalSets[i])) %in% c(8:10))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$TotalSets[i]))%in% c(1:5)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$TotalSets[i])) %ni% c(8:10)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$TotalSets[i])) %in% c(11:13))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$Group)]))))

colors13 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread %>% dplyr::select(description, color2)
colors13 <- unique(colors13)
colors13 <- colors13[order(colors13$description),]
colors13 <- colors13$color2
names(colors13) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 13.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,13,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "AA", "BB", "CC", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Spread Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors13, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread_plot ## Saved pdf 12x10 "Summary_topGO_Yazdani&Pepper&NbPCA_Up_0.7_Spread"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Spread %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbPCA_Up_0.7_Spread.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


## Down

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq=="Spread",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$TotalSets[i]))<11) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$TotalSets[i])) %in% c(1:5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$TotalSets[i])) %in% c(8:10))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$TotalSets[i]))%in% c(1:5)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$TotalSets[i])) %ni% c(8:10)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$TotalSets[i])) %in% c(11:13))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$Group)]))))

colors14 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread %>% dplyr::select(description, color2)
colors14 <- unique(colors14)
colors14 <- colors14[order(colors14$description),]
colors14 <- colors14$color2
names(colors14) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 13.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,13,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "AA", "BB", "CC", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Spread Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors14, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread_plot ## Saved pdf 12x10 "Summary_topGO_Yazdani&Pepper&NbPCA_Down_0.7_Spread

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Spread %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbPCA_Down_0.7_Spread.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


############################################################################
## Graph only Ripening-independet

## Up

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up$freq=="Ripening-independent",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$TotalSets[i]))<11)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$TotalSets[i]))%in% c(6,7)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$TotalSets[i])) %ni% c(8:10)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$TotalSets[i])) %in% c(11:13))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$Group)]))))

colors15 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind %>% dplyr::select(description, color2)
colors15 <- unique(colors15)
colors15 <- colors15[order(colors15$description),]
colors15 <- colors15$color2
names(colors15) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 13.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,13,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "AA", "BB", "CC", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Ripening-Independent Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors15, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind_plot ## Saved pdf 15X10 "Summary_topGO_Yazdani&Pepper&NbPCA_Up_0.7_Ripening_independent"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Up_Rip_ind %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbPCA_Up_0.7_Ripening_independent.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Down

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down$freq=="Ripening-independent",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$TotalSets[i]))<11)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$TotalSets[i]))%in% c(6,7)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$TotalSets[i])) %ni% c(8:10)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$TotalSets[i])) %in% c(11:13))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$Group)]))))

colors16 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind %>% dplyr::select(description, color2)
colors16 <- unique(colors16)
colors16 <- colors16[order(colors16$description),]
colors16 <- colors16$color2
names(colors16) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 13.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,13,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "AA", "BB", "CC", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Ripening-Independent Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors16, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind_plot ## Saved pdf 15x10 "Summary_topGO_Yazdani&Pepper&NbPCA_Down_0.7_Ripening_independent"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_PCA_Down_Rip_ind %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbPCA_Down_0.7_Ripening_independent.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


##################################################################################################
########################################  Distance Matrix  ###########################################
##################################################################################################
## YAZDANI VS NICOTIANA

########### Up

filelist_Yazdani_vs_Niben_DistM_Up <- lapply(list("topGO_weightFS_M82_B_Up_BP_4.0_AHRD_REVIGO",
                                                "topGO_weightFS_M82_or_Up_BP_4.0_AHRD_REVIGO",
                                                "topGO_weightFS_M82_Red_Up_BP_4.0_AHRD_REVIGO",
                                                "topGO_weightFS_MG_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                "topGO_weightFS_B_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                "topGO_weightFS_22h_Up_BP_REVIGO",
                                                "topGO_weightFS_YY_Up_BP_REVIGO",
                                                "topGO_weightFS_ZZ_Up_BP_REVIGO",
                                                "topGO_weightFS_40h_Up_BP_REVIGO",
                                                "topGO_weightFS_46h_Up_BP_REVIGO",
                                                "topGO_weightFS_56h_Up_BP_REVIGO",
                                                "topGO_weightFS_96h_Up_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Niben_DistM_Up <- do.call(rbind, filelist_Yazdani_vs_Niben_DistM_Up)

Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Stage <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Stage,
                                                        levels = c("B", "or", "Red", "ORHis", "22h", "YY", "ZZ","40h","46h","56h","96h"),
                                                        labels = c("Breaker", "Orange", "Red", "ORHis", "22h", "YY", "ZZ","40h","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Background <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Background,
                                                             levels = c("M82", "MG", "B", "Niben"),
                                                             labels = c("M82", "MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Group <- paste(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Stage, Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Group <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Group,
                                                        levels = c("Breaker_M82", "Orange_M82", "Red_M82", "ORHis_MG", "ORHis_Breaker","22h_Niben", "YY_Niben", "ZZ_Niben", "40h_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Niben_DistM_Up <- Summary_GOterms_Yazdani_vs_Niben_DistM_Up[order(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_DistM_Up)){
  Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Group[Summary_GOterms_Yazdani_vs_Niben_DistM_Up$description==Summary_GOterms_Yazdani_vs_Niben_DistM_Up$description[i]])
  Summary_GOterms_Yazdani_vs_Niben_DistM_Up$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Group[Summary_GOterms_Yazdani_vs_Niben_DistM_Up$description==Summary_GOterms_Yazdani_vs_Niben_DistM_Up$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_DistM_Up)){
  if(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$color[i] <- "#8DA0CB"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$NumberSets[i]>9){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$color[i] <- "#66C2A5"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$color[i] <- "#E78AC3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$color[i] <- "#FC8D62"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i]))>5)){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$color[i] <- "#B3B3B3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i]))>3) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i])) %in% c(6:12))){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$color[i] <- "#FFD92F"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$TotalSets[i])) %in% c(1:3))){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$color[i] <- "#E5C494"
  } else {
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Up$color[i] <- "#A6D854"}
}

## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq, levels = c("Constant", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Niben_DistM_Up$description <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$description[order(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Group, Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors17 <- Summary_GOterms_Yazdani_vs_Niben_DistM_Up %>% dplyr::select(description, color)
colors17 <- unique(colors17)
colors17 <- colors17[order(colors17$description),]
colors17 <- colors17$color
names(colors17) <- rev(unique(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$description[order(Summary_GOterms_Yazdani_vs_Niben_DistM_Up$Group, Summary_GOterms_Yazdani_vs_Niben_DistM_Up$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Niben_DistM_Up_plot <- ggplot(Summary_GOterms_Yazdani_vs_Niben_DistM_Up, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#66C2A5",
                               "Ripening"="#E78AC3",
                               "Nicotiana-independent"="#E5C494",
                               "OR-Dependent"="#FC8D62", 
                               "Ripening-independent"="#FFD92F", 
                               "Spread"="#A6D854", 
                               "Nicotiana"="#B3B3B3",
                               "Exclusive"="#8DA0CB"),
                     breaks = c("Constant",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Nicotiana Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme(axis.text.y = element_text(color = colors17, face = "bold"))
Summary_GOterms_Yazdani_vs_Niben_DistM_Up_plot ## Saved pdf 20x70 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Yazdani_vs_Niben_DistM_Up %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&NbDistM_Up_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


########### Down

filelist_Yazdani_vs_Niben_DistM_Down <- lapply(list("topGO_weightFS_M82_B_Down_BP_4.0_AHRD_REVIGO",
                                                  "topGO_weightFS_M82_or_Down_BP_4.0_AHRD_REVIGO",
                                                  "topGO_weightFS_M82_Red_Down_BP_4.0_AHRD_REVIGO",
                                                  "topGO_weightFS_MG_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                  "topGO_weightFS_B_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                  "topGO_weightFS_22h_Down_BP_REVIGO",
                                                  "topGO_weightFS_YY_Down_BP_REVIGO",
                                                  "topGO_weightFS_ZZ_Down_BP_REVIGO",
                                                  "topGO_weightFS_40h_Down_BP_REVIGO",
                                                  "topGO_weightFS_46h_Down_BP_REVIGO",
                                                  "topGO_weightFS_56h_Down_BP_REVIGO",
                                                  "topGO_weightFS_96h_Down_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Niben_DistM_Down <- do.call(rbind, filelist_Yazdani_vs_Niben_DistM_Down)

Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Stage <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Stage,
                                                          levels = c("B", "or", "Red", "ORHis", "22h", "YY", "ZZ", "40h","46h","56h","96h"),
                                                          labels = c("Breaker", "Orange", "Red", "ORHis", "22h", "YY", "ZZ","40h","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Background <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Background,
                                                               levels = c("M82", "MG", "B", "Niben"),
                                                               labels = c("M82", "MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Group <- paste(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Stage, Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Group <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Group,
                                                          levels = c("Breaker_M82", "Orange_M82", "Red_M82", "ORHis_MG", "ORHis_Breaker","22h_Niben", "YY_Niben", "ZZ_Niben", "40h_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Niben_DistM_Down <- Summary_GOterms_Yazdani_vs_Niben_DistM_Down[order(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_DistM_Down)){
  Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Group[Summary_GOterms_Yazdani_vs_Niben_DistM_Down$description==Summary_GOterms_Yazdani_vs_Niben_DistM_Down$description[i]])
  Summary_GOterms_Yazdani_vs_Niben_DistM_Down$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Group[Summary_GOterms_Yazdani_vs_Niben_DistM_Down$description==Summary_GOterms_Yazdani_vs_Niben_DistM_Down$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Niben_DistM_Down)){
  if(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$color[i] <- "#8DA0CB"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$NumberSets[i]>9){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$color[i] <- "#66C2A5"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$color[i] <- "#E78AC3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$color[i] <- "#FC8D62"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i]))>5)){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$color[i] <- "#B3B3B3"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i]))>3) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i])) %in% c(6:12))){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$color[i] <- "#FFD92F"
  } else if(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i])) %in% c(4,5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$TotalSets[i])) %in% c(1:3))){
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$color[i] <- "#E5C494"
  } else {
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Niben_DistM_Down$color[i] <- "#A6D854"}
}

## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq, levels = c("Constant", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Niben_DistM_Down$description <- factor(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$description[order(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Group, Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors18 <- Summary_GOterms_Yazdani_vs_Niben_DistM_Down %>% dplyr::select(description, color)
colors18 <- unique(colors18)
colors18 <- colors18[order(colors18$description),]
colors18 <- colors18$color
names(colors18) <- rev(unique(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$description[order(Summary_GOterms_Yazdani_vs_Niben_DistM_Down$Group, Summary_GOterms_Yazdani_vs_Niben_DistM_Down$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Niben_DistM_Down_plot <- ggplot(Summary_GOterms_Yazdani_vs_Niben_DistM_Down, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#66C2A5",
                               "Ripening"="#E78AC3",
                               "Nicotiana-independent"="#E5C494",
                               "OR-Dependent"="#FC8D62", 
                               "Ripening-independent"="#FFD92F", 
                               "Spread"="#A6D854", 
                               "Nicotiana"="#B3B3B3",
                               "Exclusive"="#8DA0CB"),
                     breaks = c("Constant",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Nicotiana Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme(axis.text.y = element_text(color = colors18, face = "bold"))
Summary_GOterms_Yazdani_vs_Niben_DistM_Down_plot ## Saved pdf 20x70 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Yazdani_vs_Niben_DistM_Down %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&NbDistM_Down_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)



##################################################################################################
## YAZDANI VS PEPPER VS NICOTIANA

########### Up

filelist_Yazdani_vs_Pepper_vs_Niben_DistM_Up <- lapply(list("topGO_weightFS_M82_B_Up_BP_4.0_AHRD_REVIGO",
                                                          "topGO_weightFS_M82_or_Up_BP_4.0_AHRD_REVIGO",
                                                          "topGO_weightFS_M82_Red_Up_BP_4.0_AHRD_REVIGO",
                                                          "topGO_weightFS_Breaker_Up_BP_REVIGO",
                                                          "topGO_weightFS_Red_Up_BP_REVIGO",
                                                          "topGO_weightFS_MG_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                          "topGO_weightFS_B_ORHis_Up_BP_4.0_AHRD_REVIGO",
                                                          "topGO_weightFS_22h_Up_BP_REVIGO",
                                                          "topGO_weightFS_YY_Up_BP_REVIGO",
                                                          "topGO_weightFS_ZZ_Up_BP_REVIGO",
                                                          "topGO_weightFS_40h_Up_BP_REVIGO",
                                                          "topGO_weightFS_46h_Up_BP_REVIGO",
                                                          "topGO_weightFS_56h_Up_BP_REVIGO",
                                                          "topGO_weightFS_96h_Up_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up <- do.call(rbind, filelist_Yazdani_vs_Pepper_vs_Niben_DistM_Up)

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Stage <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Stage,
                                                                  levels = c("B", "or", "Red", "Breaker", "Red", "ORHis", "22h", "YY", "ZZ","40h","46h","56h","96h"),
                                                                  labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "ORHis", "22h", "YY", "ZZ","40h","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Background <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Background,
                                                                       levels = c("M82", "Pepper","MG", "B", "Niben"),
                                                                       labels = c("M82", "Pepper","MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Group <- paste(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Stage, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Group <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Group,
                                                                  levels = c("Breaker_M82", "Orange_M82", "Red_M82", "Breaker_Pepper", "Red_Pepper","ORHis_MG", "ORHis_Breaker","22h_Niben", "YY_Niben", "ZZ_Niben", "40h_Niben","46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up)){
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$description[i]])
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up)){
  if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i]>11){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#370335B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "Tomato-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#bea3ce"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "Pepper-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#fccde5"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i])) %in% c(1:3)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i])) %in% c(4,5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#bc80bd"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i]))>5) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i]))<8)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#F05C3BB2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i]))>7)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#D2AF81B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i]))>5) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i])) %in% c(8:14))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#FED439B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i]))<8) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$TotalSets[i])) %in% c(1:5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#808080"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq, levels = c("Constant","Tomato-Ripening", "Pepper-Ripening", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors19 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up %>% dplyr::select(description, color)
colors19 <- unique(colors19)
colors19 <- colors19[order(colors19$description),]
colors19 <- colors19$color
names(colors19) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#370335B2",
                               "Tomato-Ripening" = "#bea3ce",#"#bebada",
                               "Pepper-Ripening" = "#fccde5",
                               "Ripening"="#bc80bd",
                               "Nicotiana-independent"="#808080",#"#D2AF81B2",
                               "OR-Dependent"="#F05C3BB2", 
                               "Ripening-independent"="#FED439B2", 
                               "Spread"="#b3de69", 
                               "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                               "Exclusive"="#197EC0B2"),
                     breaks = c("Constant",
                                "Tomato-Ripening",
                                "Pepper-Ripening",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme(axis.text.y = element_text(color = colors19, face = "bold"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_plot ## Saved pdf 20x73 "Summary_topGO_Yazdani&NbTimePoints_Up_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbDistM_Up_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

########### Down

filelist_Yazdani_vs_Pepper_vs_Niben_DistM_Down <- lapply(list("topGO_weightFS_M82_B_Down_BP_4.0_AHRD_REVIGO",
                                                            "topGO_weightFS_M82_or_Down_BP_4.0_AHRD_REVIGO",
                                                            "topGO_weightFS_M82_Red_Down_BP_4.0_AHRD_REVIGO",
                                                            "topGO_weightFS_Breaker_Down_BP_REVIGO",
                                                            "topGO_weightFS_Red_Down_BP_REVIGO",
                                                            "topGO_weightFS_MG_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                            "topGO_weightFS_B_ORHis_Down_BP_4.0_AHRD_REVIGO",
                                                            "topGO_weightFS_22h_Down_BP_REVIGO",
                                                            "topGO_weightFS_YY_Down_BP_REVIGO",
                                                            "topGO_weightFS_ZZ_Down_BP_REVIGO",
                                                            "topGO_weightFS_40h_Down_BP_REVIGO",
                                                            "topGO_weightFS_46h_Down_BP_REVIGO",
                                                            "topGO_weightFS_56h_Down_BP_REVIGO",
                                                            "topGO_weightFS_96h_Down_BP_REVIGO"), get)
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down <- do.call(rbind, filelist_Yazdani_vs_Pepper_vs_Niben_DistM_Down)

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Stage <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Stage,
                                                                    levels = c("B", "or", "Red", "Breaker", "Red", "ORHis", "22h", "YY", "ZZ","40h","46h","56h","96h"),
                                                                    labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "ORHis", "22h", "YY", "ZZ","40h","46h","56h","96h"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Background <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Background,
                                                                         levels = c("M82", "Pepper","MG", "B", "Niben"),
                                                                         labels = c("M82", "Pepper","MG", "Breaker", "Niben"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Group <- paste(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Stage, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Background, sep = "_")
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Group <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Group,
                                                                    levels = c("Breaker_M82", "Orange_M82", "Red_M82", "Breaker_Pepper", "Red_Pepper","ORHis_MG", "ORHis_Breaker","22h_Niben", "YY_Niben", "ZZ_Niben", "40h_Niben", "46h_Niben", "56h_Niben", "96h_Niben"))

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Group),]

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down)){
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i] <- list(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$description[i]])
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i] <- length(c(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Group[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$description==Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$description[i]]))
}

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down)){
  if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i]==1){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "Exclusive"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#197EC0B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i]>11){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "Constant"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#370335B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i]))<4)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "Tomato-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#bea3ce"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i]))>3) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i]))<6)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "Pepper-Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#fccde5"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i]))<6) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i])) %in% c(1:3)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i])) %in% c(4,5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "Ripening"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#bc80bd"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i]))>5) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i]))<8)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "OR-Dependent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#F05C3BB2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i]))>7)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "Nicotiana"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#D2AF81B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i]))>5) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i])) %in% c(8:14))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "Ripening-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#FED439B2"
  } else if(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$NumberSets[i]>1 & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i]))<8) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i])) %in% c(6,7)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$TotalSets[i])) %in% c(1:5))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "Nicotiana-independent"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#808080"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq[i] <- "Spread"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$color[i] <- "#b3de69"}
}


## We transform into factors freq and description variables. The last one is important for the order
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq, levels = c("Constant","Tomato-Ripening", "Pepper-Ripening", "Ripening","Nicotiana-independent", "OR-Dependent", "Ripening-independent", "Spread", "Nicotiana","Exclusive"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq)]))))

## We create a variable to assign color to the row label, in the same order as the 'description' variable
colors20 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down %>% dplyr::select(description, color)
colors20 <- unique(colors20)
colors20 <- colors20[order(colors20$description),]
colors20 <- colors20$color
names(colors20) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$Group, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq)]))


## We create the plot
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down, aes(x=Group, y=description))+
  geom_point(aes(col=freq),size=3)+
  scale_color_manual(values= c("Constant"="#370335B2",
                               "Tomato-Ripening" = "#bea3ce",#"#bebada",
                               "Pepper-Ripening" = "#fccde5",
                               "Ripening"="#bc80bd",
                               "Nicotiana-independent"="#808080",#"#D2AF81B2",
                               "OR-Dependent"="#F05C3BB2", 
                               "Ripening-independent"="#FED439B2", 
                               "Spread"="#b3de69", 
                               "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                               "Exclusive"="#197EC0B2"),
                     breaks = c("Constant",
                                "Tomato-Ripening",
                                "Pepper-Ripening",
                                "Ripening",
                                "Nicotiana-independent",
                                "OR-Dependent", 
                                "Ripening-independent", 
                                "Spread", 
                                "Nicotiana",
                                "Exclusive"))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme(axis.text.y = element_text(color = colors20, face = "bold"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_plot ## Saved pdf 20x73 "Summary_topGO_Yazdani&NbTimePoints_Down_0.7"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbDistM_Down_0.7.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

######################################################
## SUMMARIZING RESULTS BY GROUP (% AND COUNTS)

## Up
percentData_Yazdani_vs_Pepper_vs_Niben_DistM_Up <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up %>%
  group_by(Group) %>% count(freq) %>%
  mutate(ratio=scales::percent(n/sum(n)))
percentData_Yazdani_vs_Pepper_vs_Niben_DistM_Up

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_percent <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up, aes(x= Group))+
  geom_bar(aes(fill=freq), position = position_fill(reverse = TRUE))+
  geom_text(data = percentData_Yazdani_vs_Pepper_vs_Niben_DistM_Up, aes(y=n, label=ratio),
            position = position_fill(vjust = 0.5), colour = "white", fontface = "bold")+
  scale_fill_manual(values= c("Constant"="#370335B2",
                              "Tomato-Ripening" = "#bea3ce",#"#bebada",
                              "Pepper-Ripening" = "#fccde5",
                              "Ripening"="#bc80bd",
                              "Nicotiana-independent"="#808080",#"#D2AF81B2",
                              "OR-Dependent"="#F05C3BB2", 
                              "Ripening-independent"="#FED439B2", 
                              "Spread"="#b3de69", 
                              "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                              "Exclusive"="#197EC0B2"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     labels = c("0", "25", "50", "75", "100"))+
  ylab("% of GO terms subgroups")+
  theme_HPLC+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_percent 

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_abs <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up, aes(x=Group))+
  geom_bar(aes(fill=Background))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_discrete(labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "YY", "ZZ","40h", "46h", "56h", "96h"))+
  ylab("Number of GO terms")+
  scale_fill_manual(values = c("M82" = "#bea3ce",
                               "Pepper" ="#fccde5",
                               "MG"="#e45132",
                               "Breaker"="#bd2e15",
                               "Niben"="#D2AF81B2"))+
  theme_HPLC+
  theme(legend.position = "top",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_abs

figureE <- ggarrange(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_abs,
                     Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_percent,
                     ncol = 1,
                     nrow = 2)
figureE ## Saved pdf 15x10 "Summary_topGO_Yazdani&Pepper&NbDistM_Up_0.7_Percent"


## Down
percentData_Yazdani_vs_Pepper_vs_Niben_DistM_Down <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down %>%
  group_by(Group) %>% count(freq) %>%
  mutate(ratio=scales::percent(n/sum(n)))
percentData_Yazdani_vs_Pepper_vs_Niben_DistM_Down

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_percent <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down, aes(x= Group))+
  geom_bar(aes(fill=freq), position = position_fill(reverse = TRUE))+
  geom_text(data = percentData_Yazdani_vs_Pepper_vs_Niben_DistM_Down, aes(y=n, label=ratio),
            position = position_fill(vjust = 0.5), colour = "white", fontface = "bold")+
  scale_fill_manual(values= c("Constant"="#370335B2",
                              "Tomato-Ripening" = "#bea3ce",#"#bebada",
                              "Pepper-Ripening" = "#fccde5",
                              "Ripening"="#bc80bd",
                              "Nicotiana-independent"="#808080",#"#D2AF81B2",
                              "OR-Dependent"="#F05C3BB2", 
                              "Ripening-independent"="#FED439B2", 
                              "Spread"="#b3de69", 
                              "Nicotiana"="#D2AF81B2",#"#8A9197B2",
                              "Exclusive"="#197EC0B2"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     labels = c("0", "25", "50", "75", "100"))+
  ylab("% of GO terms subGroups")+
  theme_HPLC+
  theme(legend.position = "bottom",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_percent 

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_abs <- ggplot(
  Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down, aes(x=Group))+
  geom_bar(aes(fill=Background))+
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_discrete(labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "YY", "ZZ","40h", "46h", "56h", "96h"))+
  scale_fill_manual(values = c("M82" = "#bea3ce",
                               "Pepper" ="#fccde5",
                               "MG"="#e45132",
                               "Breaker"="#bd2e15",
                               "Niben"="#D2AF81B2"))+
  theme_HPLC+
  theme(legend.position = "top",
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_abs

figureF <- ggarrange(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_abs,
                     Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_percent,
                     ncol = 1,
                     nrow = 2)
figureF ## Saved pdf 15x10 "Summary_topGO_Yazdani&Pepper&NbDistM_Down_0.7_Percent"



############################################################################
## Graph only Spread

## Up

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq=="Spread",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$TotalSets[i]))<11) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$TotalSets[i])) %in% c(1:5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$TotalSets[i])) %in% c(8:10))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$TotalSets[i]))%in% c(1:5)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$TotalSets[i])) %ni% c(8:10)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$TotalSets[i])) %in% c(11:14))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$Group)]))))

colors21 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread %>% dplyr::select(description, color2)
colors21 <- unique(colors21)
colors21 <- colors21[order(colors21$description),]
colors21 <- colors21$color2
names(colors21) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 14.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,14,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "YY", "ZZ", "40h", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Spread Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors21, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread_plot ## Saved pdf 12x10 "Summary_topGO_Yazdani&Pepper&NbDistM_Up_0.7_Spread"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Spread %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbDistM_Up_0.7_Spread.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


## Down

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq=="Spread",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$TotalSets[i]))<11) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$TotalSets[i])) %in% c(1:5)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$TotalSets[i])) %in% c(8:10))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$TotalSets[i]))%in% c(1:5)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$TotalSets[i])) %ni% c(8:10)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$TotalSets[i])) %in% c(11:14))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$Group)]))))

colors22 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread %>% dplyr::select(description, color2)
colors22 <- unique(colors22)
colors22 <- colors22[order(colors22$description),]
colors22 <- colors22$color2
names(colors22) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 14.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,14,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "YY", "ZZ", "40h", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Spread Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors22, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread_plot ## Saved pdf 12x10 "Summary_topGO_Yazdani&Pepper&NbDistM_Down_0.7_Spread

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Spread %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbDistM_Down_0.7_Spread.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


############################################################################
## Graph only Ripening-independet

## Up

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up$freq=="Ripening-independent",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$TotalSets[i]))<11)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$TotalSets[i]))%in% c(6,7)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$TotalSets[i])) %ni% c(8:10)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$TotalSets[i])) %in% c(11:14))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$Group)]))))

colors23 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind %>% dplyr::select(description, color2)
colors23 <- unique(colors23)
colors23 <- colors23[order(colors23$description),]
colors23 <- colors23$color2
names(colors23) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 14.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,14,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "YY", "ZZ", "40h", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Ripening-Independent Up-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Up" = "22h",
  #                           "25h_Up" = "25h",
  #                           "28h_Up" = "28h",
  #                           "34h_Up" = "34h",
  #                           "37h_Up" = "37h",
  #                           "40h_Up" = "40h",
  #                           "46h_Up" = "46h",
  #                           "56h_Up" = "56h",
  #                           "96h_Up" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors23, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind_plot ## Saved pdf 15X10 "Summary_topGO_Yazdani&Pepper&NbDistM_Up_0.7_Ripening_independent"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Up_Rip_ind %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbDistM_Up_0.7_Ripening_independent.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## Down

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down[Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down$freq=="Ripening-independent",]

'%ni%' <- Negate('%in%')

for(i in 1:nrow(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind)){
  if(all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$TotalSets[i]))<11)){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$freq2[i] <- "Early"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$color2[i] <- "#1b9e77"
  } else if(any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$TotalSets[i]))%in% c(6,7)) & all(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$TotalSets[i])) %ni% c(8:10)) & any(as.integer(unlist(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$TotalSets[i])) %in% c(11:14))){
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$freq2[i] <- "Late"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$color2[i] <- "#d95f02"
  } else {
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$freq2[i] <- "Disperse"
    Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$color2[i] <- "#7570b3"}
}

Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$freq2 <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$freq2, levels = c("Early", "Late", "Disperse"))
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$description <- factor(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$description, levels = (rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$Group)]))))

colors24 <- Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind %>% dplyr::select(description, color2)
colors24 <- unique(colors24)
colors24 <- colors24[order(colors24$description),]
colors24 <- colors24$color2
names(colors24) <- rev(unique(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$description[order(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$freq2, Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind$Group)]))


Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind_plot <- ggplot(Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind, aes(x=as.numeric(Group), y=description))+
  annotate("rect",xmin = 0.5, xmax = 3.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bea3ce", alpha = 0.4)+
  annotate("rect",xmin = 3.5, xmax = 5.5,
           ymin = -Inf, ymax = Inf,
           fill = "#fccde5", alpha = 0.4)+
  annotate("rect",xmin = 5.5, xmax = 6.5,
           ymin = -Inf, ymax = Inf,
           fill = "#e45132", alpha = 0.4)+
  annotate("rect",xmin = 6.5, xmax = 7.5,
           ymin = -Inf, ymax = Inf,
           fill = "#bd2e15", alpha = 0.4)+
  annotate("rect",xmin = 7.5, xmax = 14.5,
           ymin = -Inf, ymax = Inf,
           fill = "#D2AF81B2", alpha = 0.4)+
  geom_point(aes(col=freq2),size=3)+
  scale_color_manual(values= c("Early"="#1b9e77",
                               "Late" = "#d95f02",
                               "Disperse" = "#7570b3"),
                     breaks = c("Early",
                                "Late",
                                "Disperse"))+
  scale_x_continuous(breaks = seq(1,14,1),
                     labels = c("Breaker", "Orange", "Red", "Breaker", "Red", "MG", "Breaker", "22h", "YY", "ZZ", "40h", "46h", "56h", "96h"),
                     expand = expansion(mult = 0, add = 0))+
  ggtitle("Biological Processes Yazdani-Pepper-Nicotiana Ripening-Independent Down-regulated 0.7")+
  #scale_x_discrete(labels= c("22h_Down" = "22h",
  #                           "25h_Down" = "25h",
  #                           "28h_Down" = "28h",
  #                           "34h_Down" = "34h",
  #                           "37h_Down" = "37h",
  #                           "40h_Down" = "40h",
  #                           "46h_Down" = "46h",
  #                           "56h_Down" = "56h",
  #                           "96h_Down" = "96h"))+
  theme_HPLC +
  theme(axis.text.y = element_text(color = colors24, face = "bold"),
        text = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, face = "bold", hjust = 1),
        axis.ticks.x = element_blank())
Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind_plot ## Saved pdf 15x10 "Summary_topGO_Yazdani&Pepper&NbDistM_Down_0.7_Ripening_independent"

write.table((Summary_GOterms_Yazdani_vs_Pepper_vs_Niben_DistM_Down_Rip_ind %>% dplyr::select(description, Group, freq)),
            "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Summary_topGO_Yazdani&Pepper&NbDistM_Down_0.7_Ripening_independent.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
