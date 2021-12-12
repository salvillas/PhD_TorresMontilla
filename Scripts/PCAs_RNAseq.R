library(DESeq2)
library(ggplot2)
library(gplots)
library(genefilter)
library(GenomicRanges)
library(plyr)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(GEOquery)
library(Biobase)
library(multiClust)
library(preprocessCore)
library(ctc)
library(gplots)
library(dendextend)
library(graphics)
library(grDevices)
library(amap)
library(survival)
library(prospectr)


## We spicify which files to read from the current wd. 
sampleFiles <- grep("counts", list.files("/Volumes/Seagate Expansion Drive/RNAseq/Cornell/04_quantification/STAR/HTSeq/Counts files/"), value = TRUE)

## We include all the features that describe the samples
sampleCondition <- read.table("/Volumes/Seagate Expansion Drive/RNAseq/Cornell/04_quantification/HISAT/Stringtie/Cornell_Data.txt", 
                              header = TRUE, sep = "\t")
sampleName <- gsub("_Niben261_STAR_HTSeq.counts", "", sampleFiles)
sampleTable <- sampleCondition
sampleTable$fileName <- sampleFiles
colnames(sampleTable)[1] <- "sampleName"
sampleTable$Construct <- factor(sampleTable$Construct, levels = c("GFP", "pcrtB"))
sampleTable$TimePoint <- as.factor(sampleTable$TimePoint) 
sampleTable$Replicate <- as.factor(sampleTable$Replicate)
str(sampleTable)

colOrder <- c("sampleName", "fileName", "Construct", "TimePoint", "Replicate")
sampleTable <- sampleTable[, colOrder]

## We build the DESeqDataSet
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "/Volumes/Seagate Expansion Drive/RNAseq/Cornell/04_quantification/STAR/HTSeq/Counts files/",
                                       design = ~Construct + TimePoint + Construct:TimePoint)
ddsHTSeq
summary(counts(ddsHTSeq))

## We pre-filter rows with 10 or less counts
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

## We create an object only with the counts
countsCornell <- counts(ddsHTSeq)
countsCornell <- as.data.frame(countsCornell)
countsCornell$Gene <- row.names(countsCornell)

###
Counts_Cornell <- countsCornell %>% select(-Gene)
Counts_Cornell <- as.matrix(Counts_Cornell)
numsamples <- c(1:48)

## Primero hay que normalizar los counts en función del tamaño de las librerías
y <- DGEList(counts= Counts_Cornell, group = numsamples)
y <- calcNormFactors(y)
z <- cpm(y, normalize.lib.size=TRUE)
scaledata <- scale(z)
scaledata <- scaledata[complete.cases(scaledata),]
colnames(scaledata)
scaledata <- as.data.frame(scaledata)



DESeq2PCA_2 <- prcomp(t(scaledata %>% select(-Gene)))
DESeq2PCA_Noblock <- as.data.frame(DESeq2PCA_2$x)
DESeq2PCA_Noblock <- cbind(DESeq2PCA_Noblock, sampleTable)
DESeq2PCA_Noblock$Construct <- factor(DESeq2PCA_Noblock$Construct, levels = c("GFP", "pcrtB"), labels = c("GFP", "(p)crtB"))
DESeq2PCA_Noblock$TimePoint <- factor(DESeq2PCA_Noblock$TimePoint, levels = c("22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h"))
DESeq2PCA_Noblock$Construct <- as.factor(DESeq2PCA_Noblock$Construct)

## Calculo el percentaje que explica cada PCA
DESeq2percentaje_Noblock <- round(DESeq2PCA_2$sdev / sum(DESeq2PCA_2$sdev) * 100, 2)
DESeq2percentaje_Noblock <- paste(colnames(DESeq2PCA_Noblock), "(", paste(as.character(DESeq2percentaje_Noblock), "%", ")", sep = ""))
DESeq2var_explained_Noblock <- round((DESeq2PCA_2$sdev^2/sum(DESeq2PCA_2$sdev^2))*100,2)
DESeq2var_explained_Noblock <- paste(colnames(DESeq2PCA_Noblock), "(", paste(as.character(DESeq2var_explained_Noblock), "%", ")", sep = ""))

DESeq2p_Noblock <- ggplot(DESeq2PCA_Noblock, aes(x= PC1, y= PC2, fill = Construct, label = row.names(DESeq2PCA_Noblock)))+
  geom_point(shape = 22, colour = "white", size = 3)+
  geom_text_repel(aes(color=TimePoint))+
  scale_color_manual(values = c("22h" ="#E64B35B2", "25h" = "#4DBBD5B2", "28h" = "#00A087B2", "34h" = "#F39B7FB2", "37h" = "#8491B4B2", "40h" = "#91D1C2B2", "46h" = "#7E6148B2", "56h" = "#B09C85B2"))+
  scale_fill_manual(values = c(GFP= "#66c2a5",
                               `(p)crtB`= "#fc8d62"))+
  xlab(DESeq2var_explained_Noblock[1])+
  ylab(DESeq2var_explained_Noblock[2])+
  theme_classic()+
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom")
DESeq2p_Noblock



HPLC_PAM_data <- read.delim("/Volumes/Seagate Expansion Drive/RNAseq/Cornell/04_quantification/HISAT/Stringtie/RNAseq_Cornell_HPLC_PAM.txt",
                            sep = "\t", header = TRUE, dec = ".")

##Adecuo el archivo al mismo formato que el de expresión
HPLC_PAM_data$Samples <- paste0(HPLC_PAM_data$Construct, HPLC_PAM_data$Time_Point, HPLC_PAM_data$Replicate)
row.names(HPLC_PAM_data) <- HPLC_PAM_data$Samples

HPLC_PAM_data <- HPLC_PAM_data %>%
  select(Phytoene, Total_Carotenoids, PAM)

HPLC_PAM_data <- t(HPLC_PAM_data)
HPLC_PAM_data <- as.data.frame(HPLC_PAM_data)

## Escalo igual que la expresión
HPLC_PAM_data_scaled <- scale(HPLC_PAM_data)

################################################
#### INCLUDE HPLC BLOCK SCALED

scaledata_block <- blockScale(t(scaledata), type = "hard")
HPLC_PAM_data_block_scaled <- blockScale(t(HPLC_PAM_data_scaled), type = "hard")

DESeq2_HPLC_PAM_block <- rbind(t(scaledata_block$Xscaled), t(HPLC_PAM_data_block_scaled$Xscaled))
DESeq2_HPLC_PAM_block <- as.data.frame(DESeq2_HPLC_PAM_block)

DESeq2PCA <- prcomp(t(DESeq2_HPLC_PAM_block))
DESeq2PCA_block <- as.data.frame(DESeq2PCA$x)
DESeq2PCA_block <- cbind(DESeq2PCA_block, sampleTable)
DESeq2PCA_block$Construct <- factor(DESeq2PCA_block$Construct, levels = c("GFP", "pcrtB"), labels = c("GFP", "(p)crtB"))
DESeq2PCA_block$TimePoint <- factor(DESeq2PCA_block$TimePoint, levels = c("22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h"))


## Calculo el percentaje que explica cada PCA
DESeq2percentaje <- round(DESeq2PCA$sdev / sum(DESeq2PCA$sdev) * 100, 2)
DESeq2percentaje <- paste(colnames(DESeq2PCA_block), "(", paste(as.character(DESeq2percentaje), "%", ")", sep = ""))
DESeq2var_explained <- round((DESeq2PCA$sdev^2/sum(DESeq2PCA$sdev^2))*100,2)
DESeq2var_explained <- paste(colnames(DESeq2PCA_block), "(", paste(as.character(DESeq2var_explained), "%", ")", sep = ""))

DESeq2p_block <- ggplot(DESeq2PCA_block, aes(x= PC1, y= PC2, fill = Construct, label = row.names(DESeq2PCA_Noblock)))+
  geom_point(shape = 22, colour = "white", size = 3)+
  geom_text_repel(aes(color=TimePoint))+
  scale_color_manual(values = c("22h" ="#E64B35B2", "25h" = "#4DBBD5B2", "28h" = "#00A087B2", "34h" = "#F39B7FB2", "37h" = "#8491B4B2", "40h" = "#91D1C2B2", "46h" = "#7E6148B2", "56h" = "#B09C85B2"))+
  scale_fill_manual(values = c(GFP= "#66c2a5",
                               `(p)crtB`= "#fc8d62"))+
  xlab(DESeq2var_explained_Noblock[1])+
  ylab(DESeq2var_explained_Noblock[2])+
  theme_classic()+
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom")
DESeq2p_block


#########################################################
#### INCLUDE 96h

## We spicify which files to read from the current wd. 
sampleFiles96h <- grep("4d",grep("counts", list.files("/Volumes/Seagate Expansion Drive/RNAseq/Briardo/04_quantification/STAR/HTSeq/Counts files/"), value = TRUE), value = TRUE)

## We include all the features that describe the samples
sampleCondition2 <- read.table("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/Briardo/04_quantification/STAR/HTSeq/Briardo_data.txt", 
                               header = TRUE, sep = "\t")
colnames(sampleCondition2)[1] <- "sampleName"
sampleCondition2$fileName <- paste0(sampleCondition2$sampleName, "_Niben261_STAR_HTSeq.counts")
sampleTable2 <- sampleCondition2
sampleTable2$Construct <- factor(sampleTable2$Construct, levels = c("GFP", "pcrtB"))
sampleTable2$Replicate <- as.factor(sampleTable2$Replicate)
str(sampleTable2)


colOrder <- c("sampleName", "fileName", "Construct", "TimePoint", "Replicate")
sampleTable2 <- sampleTable2[, colOrder]

## We build the DESeqDataSet
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable2,
                                       directory = "/Volumes/Seagate Expansion Drive/RNAseq/Briardo/04_quantification/STAR/HTSeq/Counts files/",
                                       design = ~Construct)
ddsHTSeq
summary(counts(ddsHTSeq))

## We pre-filter rows with 10 or less counts
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

## We create an object only with the counts
countsBriardo <- counts(ddsHTSeq)
countsBriardo <- as.data.frame(countsBriardo)
colnames(countsBriardo) <- c("GFP_96h_1", "GFP_96h_2", "GFP_96h_3", "pcrtB_96h_1", "pcrtB_96h_2", "pcrtB_96h_3")
countsBriardo$Gene <- row.names(countsBriardo)

###
Counts_Briardo <- countsBriardo %>% select(-Gene)
Counts_Briardo <- as.matrix(Counts_Briardo)
numsamples <- c(1:6)

## Primero hay que normalizar los counts en función del tamaño de las librerías
y2 <- DGEList(counts= Counts_Briardo, group = numsamples)
y2 <- calcNormFactors(y2)
z2 <- cpm(y2, normalize.lib.size=TRUE)
scaledata2 <- scale(z2)
scaledata2 <- scaledata2[complete.cases(scaledata2),]
colnames(scaledata2)
scaledata2 <- as.data.frame(scaledata2)

scaledata_Cornell_Briardo <- merge(scaledata, scaledata2, by= 0)
rownames(scaledata_Cornell_Briardo) <- scaledata_Cornell_Briardo$Row.names
scaledata_Cornell_Briardo$Row.names <- NULL

sampleCondition2 <- sampleCondition2 %>% select(colnames(sampleTable))
sampleTable_Cornell_Briardo <- rbind(sampleTable, sampleCondition2)

DESeq2PCA_Cornell_Briardo <- prcomp(t(scaledata_Cornell_Briardo))
dfCB <- as.data.frame(DESeq2PCA_Cornell_Briardo$x)

dfCB <- cbind(dfCB, sampleTable_Cornell_Briardo)
dfCB$Construct <- factor(dfCB$Construct, levels = c("GFP", "pcrtB"), labels = c("GFP", "(p)crtB"))
dfCB$TimePoint <- factor(dfCB$TimePoint, levels = c("22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h"))

## Calculo el percentaje que explica cada PCA
DESeq2percentaje_Noblock_Cornell_Briardo <- round(DESeq2PCA_Cornell_Briardo$sdev / sum(DESeq2PCA_Cornell_Briardo$sdev) * 100, 2)
DESeq2percentaje_Noblock_Cornell_Briardo <- paste(colnames(dfCB), "(", paste(as.character(DESeq2percentaje_Noblock_Cornell_Briardo), "%", ")", sep = ""))
DESeq2var_explained_Noblock_Cornell_Briardo <- round((DESeq2PCA_Cornell_Briardo$sdev^2/sum(DESeq2PCA_2$sdev^2))*100,2)
DESeq2var_explained_Noblock_Cornell_Briardo <- paste(colnames(dfCB), "(", paste(as.character(DESeq2var_explained_Noblock_Cornell_Briardo), "%", ")", sep = ""))

DESeq2p_Noblock_Cornell_Briardo <- ggplot(dfCB, aes(x= PC1, y= PC2, fill = Construct, label = row.names(dfCB)))+
  geom_point(shape = 22, colour = "white", size = 3)+
  geom_text_repel(aes(color=TimePoint))+
  scale_color_manual(values = c("22h" ="#E64B35B2", "25h" = "#4DBBD5B2", "28h" = "#00A087B2", "34h" = "#F39B7FB2", "37h" = "#8491B4B2", "40h" = "#91D1C2B2", "46h" = "#7E6148B2", "56h" = "#B09C85B2", "96h"= "red"))+
  scale_fill_manual(values = c(GFP= "#66c2a5",
                               `(p)crtB`= "#fc8d62"))+
  xlab(DESeq2var_explained_Noblock[1])+
  ylab(DESeq2var_explained_Noblock[2])+
  theme_classic()+
  theme(legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom")
DESeq2p_Noblock_Cornell_Briardo



