library(DESeq2)
library(ggplot2)
library(gplots)
library(genefilter)
library(GenomicRanges)
library(tidyverse)
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
library(VennDiagram)
library(UpSetR)
library(ComplexHeatmap)
library(topGO)
library(seqinr)

#####################################################################################
########################### NICOTIANA CORNELL #######################################
#####################################################################################

## First, we upload the DEGs files (mapped to Niben261 and analyzed with DESeq2)

setwd("E:/RNAseq/Cornell/06_GeneOntology/Cornell")

## We use the files created in the Cornell_DESeq2.R Script

Cornell_22h_Up <- as.data.frame(resSig_22h_2[resSig_22h_2$log2FoldChange>0,])
Cornell_22h_Up$Niben <- paste0(rownames(Cornell_22h_Up),".1")
Cornell_22h_Down <- as.data.frame(resSig_22h_2[resSig_22h_2$log2FoldChange<0,])
Cornell_22h_Down$Niben <- paste0(rownames(Cornell_22h_Down),".1")

Cornell_25h_Up <- as.data.frame(resSig_25h_2[resSig_25h_2$log2FoldChange>0,])
Cornell_25h_Up$Niben <- paste0(rownames(Cornell_25h_Up),".1")
Cornell_25h_Down <- as.data.frame(resSig_25h_2[resSig_25h_2$log2FoldChange<0,])
Cornell_25h_Down$Niben <- paste0(rownames(Cornell_25h_Down),".1")

Cornell_28h_Up <- as.data.frame(resSig_28h_2[resSig_28h_2$log2FoldChange>0,])
Cornell_28h_Up$Niben <- paste0(rownames(Cornell_28h_Up),".1")
Cornell_28h_Down <- as.data.frame(resSig_28h_2[resSig_28h_2$log2FoldChange<0,])
Cornell_28h_Down$Niben <- paste0(rownames(Cornell_28h_Down),".1")

Cornell_34h_Up <- as.data.frame(resSig_34h_2[resSig_34h_2$log2FoldChange>0,])
Cornell_34h_Up$Niben <- paste0(rownames(Cornell_34h_Up),".1")
Cornell_34h_Down <- as.data.frame(resSig_34h_2[resSig_34h_2$log2FoldChange<0,])
Cornell_34h_Down$Niben <- paste0(rownames(Cornell_34h_Down),".1")

Cornell_37h_Up <- as.data.frame(resSig_37h_2[resSig_37h_2$log2FoldChange>0,])
Cornell_37h_Up$Niben <- paste0(rownames(Cornell_37h_Up),".1")
Cornell_37h_Down <- as.data.frame(resSig_37h_2[resSig_37h_2$log2FoldChange<0,])
Cornell_37h_Down$Niben <- paste0(rownames(Cornell_37h_Down),".1")

Cornell_40h_Up <- as.data.frame(resSig_40h_2[resSig_40h_2$log2FoldChange>0,])
Cornell_40h_Up$Niben <- paste0(rownames(Cornell_40h_Up),".1")
Cornell_40h_Down <- as.data.frame(resSig_40h_2[resSig_40h_2$log2FoldChange<0,])
Cornell_40h_Down$Niben <- paste0(rownames(Cornell_40h_Down),".1")

Cornell_46h_Up <- as.data.frame(resSig_46h_2[resSig_46h_2$log2FoldChange>0,])
Cornell_46h_Up$Niben <- paste0(rownames(Cornell_46h_Up),".1")
Cornell_46h_Down <- as.data.frame(resSig_46h_2[resSig_46h_2$log2FoldChange<0,])
Cornell_46h_Down$Niben <- paste0(rownames(Cornell_46h_Down),".1")

Cornell_56h_Up <- as.data.frame(resSig_56h_2[resSig_56h_2$log2FoldChange>0,])
Cornell_56h_Up$Niben <- paste0(rownames(Cornell_56h_Up),".1")
Cornell_56h_Down <- as.data.frame(resSig_56h_2[resSig_56h_2$log2FoldChange<0,])
Cornell_56h_Down$Niben <- paste0(rownames(Cornell_56h_Down),".1")

Cornell_AA_Up <- as.data.frame(resSig_AA_2[resSig_AA_2$log2FoldChange>0,])
Cornell_AA_Up$Niben <- paste0(rownames(Cornell_AA_Up),".1")
Cornell_AA_Down <- as.data.frame(resSig_AA_2[resSig_AA_2$log2FoldChange<0,])
Cornell_AA_Down$Niben <- paste0(rownames(Cornell_AA_Down),".1")

Cornell_BB_Up <- as.data.frame(resSig_BB_2[resSig_BB_2$log2FoldChange>0,])
Cornell_BB_Up$Niben <- paste0(rownames(Cornell_BB_Up),".1")
Cornell_BB_Down <- as.data.frame(resSig_BB_2[resSig_BB_2$log2FoldChange<0,])
Cornell_BB_Down$Niben <- paste0(rownames(Cornell_BB_Down),".1")

Cornell_CC_Up <- as.data.frame(resSig_CC_2[resSig_CC_2$log2FoldChange>0,])
Cornell_CC_Up$Niben <- paste0(rownames(Cornell_CC_Up),".1")
Cornell_CC_Down <- as.data.frame(resSig_CC_2[resSig_CC_2$log2FoldChange<0,])
Cornell_CC_Down$Niben <- paste0(rownames(Cornell_CC_Down),".1")

Cornell_DD_Up <- as.data.frame(resSig_DD_2[resSig_DD_2$log2FoldChange>0,])
Cornell_DD_Up$Niben <- paste0(rownames(Cornell_DD_Up),".1")
Cornell_DD_Down <- as.data.frame(resSig_DD_2[resSig_DD_2$log2FoldChange<0,])
Cornell_DD_Down$Niben <- paste0(rownames(Cornell_DD_Down),".1")

Cornell_EE_Up <- as.data.frame(resSig_EE_2[resSig_EE_2$log2FoldChange>0,])
Cornell_EE_Up$Niben <- paste0(rownames(Cornell_EE_Up),".1")
Cornell_EE_Down <- as.data.frame(resSig_EE_2[resSig_EE_2$log2FoldChange<0,])
Cornell_EE_Down$Niben <- paste0(rownames(Cornell_EE_Down),".1")

Cornell_YY_Up <- as.data.frame(resSig_YY_2[resSig_YY_2$log2FoldChange>0,])
Cornell_YY_Up$Niben <- paste0(rownames(Cornell_YY_Up),".1")
Cornell_YY_Down <- as.data.frame(resSig_YY_2[resSig_YY_2$log2FoldChange<0,])
Cornell_YY_Down$Niben <- paste0(rownames(Cornell_YY_Down),".1")

Cornell_ZZ_Up <- as.data.frame(resSig_ZZ_2[resSig_ZZ_2$log2FoldChange>0,])
Cornell_ZZ_Up$Niben <- paste0(rownames(Cornell_ZZ_Up),".1")
Cornell_ZZ_Down <- as.data.frame(resSig_ZZ_2[resSig_ZZ_2$log2FoldChange<0,])
Cornell_ZZ_Down$Niben <- paste0(rownames(Cornell_ZZ_Down),".1")

## Nicotiana 261 GO Annotation

geneID2GO_Niben261 <- readMappings("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/Niben261_genome_complete01.gene_models.go_annotations.txt")

geneNames_Niben261 <- names(geneID2GO_Niben261)

##################################################################################################
################################  TIME POINTS   ##################################################

#############################  BIOLOGICAL PROCESS  ###############################################
#################################################################################################
## 22h Up-regulated

geneList_22h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_22h_Up$Niben))
names(geneList_22h_Up) <- geneNames_Niben261

GOData_22h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_22h_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_22h_Up_BP <- runTest(GOData_22h_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_22h_Up_BP <- runTest(GOData_22h_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_22h_Up_BP <- runTest(GOData_22h_Up_BP, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_22h_Up_BP <- runTest(GOData_22h_Up_BP, algorithm = "weight01",
                                      statistic = "fisher")

allGO_22h_Up_BP <- usedGO(GOData_22h_Up_BP)
allRes_22h_Up_BP <- GenTable(GOData_22h_Up_BP, 
                            classicFisher = resultsFis_Class_22h_Up_BP,
                            weightFisher = resultsFis_Weight_22h_Up_BP,
                            #classicKS = resultsKS_Class_22h_Up_BP,
                            #weightKS = resultsKS_Weight_22h_Up_BP,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_22h_Up_BP))

pVals_Weight_22h_Up_BP <- score(resultsFis_Weight_22h_Up_BP)[score(resultsFis_Weight_22h_Up_BP)<=0.05]

GOData_22h_Up_BP_DE <- termStat(object = GOData_22h_Up_BP, whichGO = names(pVals_Weight_22h_Up_BP))
GOData_22h_Up_BP_DE$DEG <- GOData_22h_Up_BP_DE$Significant
GOData_22h_Up_BP_DE$Term <- Term(rownames(GOData_22h_Up_BP_DE))

ggplot(GOData_22h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_22h_Up_BP_DE,
            "topGO_weightFS_22h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 22h Down-regulated

geneList_22h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_22h_Down$Niben))
names(geneList_22h_Down) <- geneNames_Niben261

GOData_22h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_22h_Down, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_22h_Down_BP <- runTest(GOData_22h_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_22h_Down_BP <- runTest(GOData_22h_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_22h_Down_BP <- runTest(GOData_22h_Down_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_22h_Down_BP <- runTest(GOData_22h_Down_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_22h_Down_BP <- usedGO(GOData_22h_Down_BP)
allRes_22h_Down_BP <- GenTable(GOData_22h_Down_BP, 
                             classicFisher = resultsFis_Class_22h_Down_BP,
                             weightFisher = resultsFis_Weight_22h_Down_BP,
                             #classicKS = resultsKS_Class_22h_Down_BP,
                             #weightKS = resultsKS_Weight_22h_Down_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_22h_Down_BP))

pVals_Weight_22h_Down_BP <- score(resultsFis_Weight_22h_Down_BP)[score(resultsFis_Weight_22h_Down_BP)<=0.05]

GOData_22h_Down_BP_DE <- termStat(object = GOData_22h_Down_BP, whichGO = names(pVals_Weight_22h_Down_BP))
GOData_22h_Down_BP_DE$DEG <- GOData_22h_Down_BP_DE$Significant
GOData_22h_Down_BP_DE$Term <- Term(rownames(GOData_22h_Down_BP_DE))

ggplot(GOData_22h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_22h_Down_BP_DE,
            "topGO_weightFS_22h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 25h Up-regulated

geneList_25h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_25h_Up$Niben))
names(geneList_25h_Up) <- geneNames_Niben261

GOData_25h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_25h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_25h_Up_BP <- runTest(GOData_25h_Up_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_25h_Up_BP <- runTest(GOData_25h_Up_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_25h_Up_BP <- runTest(GOData_25h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_25h_Up_BP <- runTest(GOData_25h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_25h_Up_BP <- usedGO(GOData_25h_Up_BP)
allRes_25h_Up_BP <- GenTable(GOData_25h_Up_BP, 
                             classicFisher = resultsFis_Class_25h_Up_BP,
                             weightFisher = resultsFis_Weight_25h_Up_BP,
                             #classicKS = resultsKS_Class_25h_Up_BP,
                             #weightKS = resultsKS_Weight_25h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_25h_Up_BP))

pVals_Weight_25h_Up_BP <- score(resultsFis_Weight_25h_Up_BP)[score(resultsFis_Weight_25h_Up_BP)<=0.05]

GOData_25h_Up_BP_DE <- termStat(object = GOData_25h_Up_BP, whichGO = names(pVals_Weight_25h_Up_BP))
GOData_25h_Up_BP_DE$DEG <- GOData_25h_Up_BP_DE$Significant
GOData_25h_Up_BP_DE$Term <- Term(rownames(GOData_25h_Up_BP_DE))

ggplot(GOData_25h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_25h_Up_BP_DE,
            "topGO_weightFS_25h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 25h Down-regulated

geneList_25h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_25h_Down$Niben))
names(geneList_25h_Down) <- geneNames_Niben261

GOData_25h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_25h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_25h_Down_BP <- runTest(GOData_25h_Down_BP, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_25h_Down_BP <- runTest(GOData_25h_Down_BP, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_25h_Down_BP <- runTest(GOData_25h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_25h_Down_BP <- runTest(GOData_25h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_25h_Down_BP <- usedGO(GOData_25h_Down_BP)
allRes_25h_Down_BP <- GenTable(GOData_25h_Down_BP, 
                               classicFisher = resultsFis_Class_25h_Down_BP,
                               weightFisher = resultsFis_Weight_25h_Down_BP,
                               #classicKS = resultsKS_Class_25h_Down_BP,
                               #weightKS = resultsKS_Weight_25h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_25h_Down_BP))

pVals_Weight_25h_Down_BP <- score(resultsFis_Weight_25h_Down_BP)[score(resultsFis_Weight_25h_Down_BP)<=0.05]

GOData_25h_Down_BP_DE <- termStat(object = GOData_25h_Down_BP, whichGO = names(pVals_Weight_25h_Down_BP))
GOData_25h_Down_BP_DE$DEG <- GOData_25h_Down_BP_DE$Significant
GOData_25h_Down_BP_DE$Term <- Term(rownames(GOData_25h_Down_BP_DE))

ggplot(GOData_25h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_25h_Down_BP_DE,
            "topGO_weightFS_25h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 28h Up-regulated

geneList_28h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_28h_Up$Niben))
names(geneList_28h_Up) <- geneNames_Niben261

GOData_28h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_28h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_28h_Up_BP <- runTest(GOData_28h_Up_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_28h_Up_BP <- runTest(GOData_28h_Up_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_28h_Up_BP <- runTest(GOData_28h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_28h_Up_BP <- runTest(GOData_28h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_28h_Up_BP <- usedGO(GOData_28h_Up_BP)
allRes_28h_Up_BP <- GenTable(GOData_28h_Up_BP, 
                             classicFisher = resultsFis_Class_28h_Up_BP,
                             weightFisher = resultsFis_Weight_28h_Up_BP,
                             #classicKS = resultsKS_Class_28h_Up_BP,
                             #weightKS = resultsKS_Weight_28h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_28h_Up_BP))

pVals_Weight_28h_Up_BP <- score(resultsFis_Weight_28h_Up_BP)[score(resultsFis_Weight_28h_Up_BP)<=0.05]

GOData_28h_Up_BP_DE <- termStat(object = GOData_28h_Up_BP, whichGO = names(pVals_Weight_28h_Up_BP))
GOData_28h_Up_BP_DE$DEG <- GOData_28h_Up_BP_DE$Significant
GOData_28h_Up_BP_DE$Term <- Term(rownames(GOData_28h_Up_BP_DE))

ggplot(GOData_28h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_28h_Up_BP_DE,
            "topGO_weightFS_28h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 28h Down-regulated

geneList_28h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_28h_Down$Niben))
names(geneList_28h_Down) <- geneNames_Niben261

GOData_28h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_28h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_28h_Down_BP <- runTest(GOData_28h_Down_BP, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_28h_Down_BP <- runTest(GOData_28h_Down_BP, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_28h_Down_BP <- runTest(GOData_28h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_28h_Down_BP <- runTest(GOData_28h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_28h_Down_BP <- usedGO(GOData_28h_Down_BP)
allRes_28h_Down_BP <- GenTable(GOData_28h_Down_BP, 
                               classicFisher = resultsFis_Class_28h_Down_BP,
                               weightFisher = resultsFis_Weight_28h_Down_BP,
                               #classicKS = resultsKS_Class_28h_Down_BP,
                               #weightKS = resultsKS_Weight_28h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_28h_Down_BP))

pVals_Weight_28h_Down_BP <- score(resultsFis_Weight_28h_Down_BP)[score(resultsFis_Weight_28h_Down_BP)<=0.05]

GOData_28h_Down_BP_DE <- termStat(object = GOData_28h_Down_BP, whichGO = names(pVals_Weight_28h_Down_BP))
GOData_28h_Down_BP_DE$DEG <- GOData_28h_Down_BP_DE$Significant
GOData_28h_Down_BP_DE$Term <- Term(rownames(GOData_28h_Down_BP_DE))

ggplot(GOData_28h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_28h_Down_BP_DE,
            "topGO_weightFS_28h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 34h Up-regulated

geneList_34h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_34h_Up$Niben))
names(geneList_34h_Up) <- geneNames_Niben261

GOData_34h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_34h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_34h_Up_BP <- runTest(GOData_34h_Up_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_34h_Up_BP <- runTest(GOData_34h_Up_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_34h_Up_BP <- runTest(GOData_34h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_34h_Up_BP <- runTest(GOData_34h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_34h_Up_BP <- usedGO(GOData_34h_Up_BP)
allRes_34h_Up_BP <- GenTable(GOData_34h_Up_BP, 
                             classicFisher = resultsFis_Class_34h_Up_BP,
                             weightFisher = resultsFis_Weight_34h_Up_BP,
                             #classicKS = resultsKS_Class_34h_Up_BP,
                             #weightKS = resultsKS_Weight_34h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_34h_Up_BP))

pVals_Weight_34h_Up_BP <- score(resultsFis_Weight_34h_Up_BP)[score(resultsFis_Weight_34h_Up_BP)<=0.05]

GOData_34h_Up_BP_DE <- termStat(object = GOData_34h_Up_BP, whichGO = names(pVals_Weight_34h_Up_BP))
GOData_34h_Up_BP_DE$DEG <- GOData_34h_Up_BP_DE$Significant
GOData_34h_Up_BP_DE$Term <- Term(rownames(GOData_34h_Up_BP_DE))

ggplot(GOData_34h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_34h_Up_BP_DE,
            "topGO_weightFS_34h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 34h Down-regulated

geneList_34h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_34h_Down$Niben))
names(geneList_34h_Down) <- geneNames_Niben261

GOData_34h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_34h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_34h_Down_BP <- runTest(GOData_34h_Down_BP, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_34h_Down_BP <- runTest(GOData_34h_Down_BP, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_34h_Down_BP <- runTest(GOData_34h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_34h_Down_BP <- runTest(GOData_34h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_34h_Down_BP <- usedGO(GOData_34h_Down_BP)
allRes_34h_Down_BP <- GenTable(GOData_34h_Down_BP, 
                               classicFisher = resultsFis_Class_34h_Down_BP,
                               weightFisher = resultsFis_Weight_34h_Down_BP,
                               #classicKS = resultsKS_Class_34h_Down_BP,
                               #weightKS = resultsKS_Weight_34h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_34h_Down_BP))

pVals_Weight_34h_Down_BP <- score(resultsFis_Weight_34h_Down_BP)[score(resultsFis_Weight_34h_Down_BP)<=0.05]

GOData_34h_Down_BP_DE <- termStat(object = GOData_34h_Down_BP, whichGO = names(pVals_Weight_34h_Down_BP))
GOData_34h_Down_BP_DE$DEG <- GOData_34h_Down_BP_DE$Significant
GOData_34h_Down_BP_DE$Term <- Term(rownames(GOData_34h_Down_BP_DE))

ggplot(GOData_34h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_34h_Down_BP_DE,
            "topGO_weightFS_34h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 37h Up-regulated

geneList_37h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_37h_Up$Niben))
names(geneList_37h_Up) <- geneNames_Niben261

GOData_37h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_37h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_37h_Up_BP <- runTest(GOData_37h_Up_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_37h_Up_BP <- runTest(GOData_37h_Up_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_37h_Up_BP <- runTest(GOData_37h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_37h_Up_BP <- runTest(GOData_37h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_37h_Up_BP <- usedGO(GOData_37h_Up_BP)
allRes_37h_Up_BP <- GenTable(GOData_37h_Up_BP, 
                             classicFisher = resultsFis_Class_37h_Up_BP,
                             weightFisher = resultsFis_Weight_37h_Up_BP,
                             #classicKS = resultsKS_Class_37h_Up_BP,
                             #weightKS = resultsKS_Weight_37h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_37h_Up_BP))

pVals_Weight_37h_Up_BP <- score(resultsFis_Weight_37h_Up_BP)[score(resultsFis_Weight_37h_Up_BP)<=0.05]

GOData_37h_Up_BP_DE <- termStat(object = GOData_37h_Up_BP, whichGO = names(pVals_Weight_37h_Up_BP))
GOData_37h_Up_BP_DE$DEG <- GOData_37h_Up_BP_DE$Significant
GOData_37h_Up_BP_DE$Term <- Term(rownames(GOData_37h_Up_BP_DE))

ggplot(GOData_37h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_37h_Up_BP_DE,
            "topGO_weightFS_37h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 37h Down-regulated

geneList_37h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_37h_Down$Niben))
names(geneList_37h_Down) <- geneNames_Niben261

GOData_37h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_37h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_37h_Down_BP <- runTest(GOData_37h_Down_BP, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_37h_Down_BP <- runTest(GOData_37h_Down_BP, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_37h_Down_BP <- runTest(GOData_37h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_37h_Down_BP <- runTest(GOData_37h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_37h_Down_BP <- usedGO(GOData_37h_Down_BP)
allRes_37h_Down_BP <- GenTable(GOData_37h_Down_BP, 
                               classicFisher = resultsFis_Class_37h_Down_BP,
                               weightFisher = resultsFis_Weight_37h_Down_BP,
                               #classicKS = resultsKS_Class_37h_Down_BP,
                               #weightKS = resultsKS_Weight_37h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_37h_Down_BP))

pVals_Weight_37h_Down_BP <- score(resultsFis_Weight_37h_Down_BP)[score(resultsFis_Weight_37h_Down_BP)<=0.05]

GOData_37h_Down_BP_DE <- termStat(object = GOData_37h_Down_BP, whichGO = names(pVals_Weight_37h_Down_BP))
GOData_37h_Down_BP_DE$DEG <- GOData_37h_Down_BP_DE$Significant
GOData_37h_Down_BP_DE$Term <- Term(rownames(GOData_37h_Down_BP_DE))

ggplot(GOData_37h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_37h_Down_BP_DE,
            "topGO_weightFS_37h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 40h Up-regulated

geneList_40h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_40h_Up$Niben))
names(geneList_40h_Up) <- geneNames_Niben261

GOData_40h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_40h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_40h_Up_BP <- runTest(GOData_40h_Up_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_40h_Up_BP <- runTest(GOData_40h_Up_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_40h_Up_BP <- runTest(GOData_40h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_40h_Up_BP <- runTest(GOData_40h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_40h_Up_BP <- usedGO(GOData_40h_Up_BP)
allRes_40h_Up_BP <- GenTable(GOData_40h_Up_BP, 
                             classicFisher = resultsFis_Class_40h_Up_BP,
                             weightFisher = resultsFis_Weight_40h_Up_BP,
                             #classicKS = resultsKS_Class_40h_Up_BP,
                             #weightKS = resultsKS_Weight_40h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_40h_Up_BP))

pVals_Weight_40h_Up_BP <- score(resultsFis_Weight_40h_Up_BP)[score(resultsFis_Weight_40h_Up_BP)<=0.05]

GOData_40h_Up_BP_DE <- termStat(object = GOData_40h_Up_BP, whichGO = names(pVals_Weight_40h_Up_BP))
GOData_40h_Up_BP_DE$DEG <- GOData_40h_Up_BP_DE$Significant
GOData_40h_Up_BP_DE$Term <- Term(rownames(GOData_40h_Up_BP_DE))

ggplot(GOData_40h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_40h_Up_BP_DE,
            "topGO_weightFS_40h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 40h Down-regulated

geneList_40h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_40h_Down$Niben))
names(geneList_40h_Down) <- geneNames_Niben261

GOData_40h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_40h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_40h_Down_BP <- runTest(GOData_40h_Down_BP, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_40h_Down_BP <- runTest(GOData_40h_Down_BP, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_40h_Down_BP <- runTest(GOData_40h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_40h_Down_BP <- runTest(GOData_40h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_40h_Down_BP <- usedGO(GOData_40h_Down_BP)
allRes_40h_Down_BP <- GenTable(GOData_40h_Down_BP, 
                               classicFisher = resultsFis_Class_40h_Down_BP,
                               weightFisher = resultsFis_Weight_40h_Down_BP,
                               #classicKS = resultsKS_Class_40h_Down_BP,
                               #weightKS = resultsKS_Weight_40h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_40h_Down_BP))

pVals_Weight_40h_Down_BP <- score(resultsFis_Weight_40h_Down_BP)[score(resultsFis_Weight_40h_Down_BP)<=0.05]

GOData_40h_Down_BP_DE <- termStat(object = GOData_40h_Down_BP, whichGO = names(pVals_Weight_40h_Down_BP))
GOData_40h_Down_BP_DE$DEG <- GOData_40h_Down_BP_DE$Significant
GOData_40h_Down_BP_DE$Term <- Term(rownames(GOData_40h_Down_BP_DE))

ggplot(GOData_40h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_40h_Down_BP_DE,
            "topGO_weightFS_40h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 46h Up-regulated

geneList_46h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_46h_Up$Niben))
names(geneList_46h_Up) <- geneNames_Niben261

GOData_46h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_46h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_46h_Up_BP <- runTest(GOData_46h_Up_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_46h_Up_BP <- runTest(GOData_46h_Up_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_46h_Up_BP <- runTest(GOData_46h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_46h_Up_BP <- runTest(GOData_46h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_46h_Up_BP <- usedGO(GOData_46h_Up_BP)
allRes_46h_Up_BP <- GenTable(GOData_46h_Up_BP, 
                             classicFisher = resultsFis_Class_46h_Up_BP,
                             weightFisher = resultsFis_Weight_46h_Up_BP,
                             #classicKS = resultsKS_Class_46h_Up_BP,
                             #weightKS = resultsKS_Weight_46h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_46h_Up_BP))

pVals_Weight_46h_Up_BP <- score(resultsFis_Weight_46h_Up_BP)[score(resultsFis_Weight_46h_Up_BP)<=0.05]

GOData_46h_Up_BP_DE <- termStat(object = GOData_46h_Up_BP, whichGO = names(pVals_Weight_46h_Up_BP))
GOData_46h_Up_BP_DE$DEG <- GOData_46h_Up_BP_DE$Significant
GOData_46h_Up_BP_DE$Term <- Term(rownames(GOData_46h_Up_BP_DE))

ggplot(GOData_46h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_46h_Up_BP_DE,
            "topGO_weightFS_46h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 46h Down-regulated

geneList_46h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_46h_Down$Niben))
names(geneList_46h_Down) <- geneNames_Niben261

GOData_46h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_46h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_46h_Down_BP <- runTest(GOData_46h_Down_BP, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_46h_Down_BP <- runTest(GOData_46h_Down_BP, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_46h_Down_BP <- runTest(GOData_46h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_46h_Down_BP <- runTest(GOData_46h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_46h_Down_BP <- usedGO(GOData_46h_Down_BP)
allRes_46h_Down_BP <- GenTable(GOData_46h_Down_BP, 
                               classicFisher = resultsFis_Class_46h_Down_BP,
                               weightFisher = resultsFis_Weight_46h_Down_BP,
                               #classicKS = resultsKS_Class_46h_Down_BP,
                               #weightKS = resultsKS_Weight_46h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_46h_Down_BP))

pVals_Weight_46h_Down_BP <- score(resultsFis_Weight_46h_Down_BP)[score(resultsFis_Weight_46h_Down_BP)<=0.05]

GOData_46h_Down_BP_DE <- termStat(object = GOData_46h_Down_BP, whichGO = names(pVals_Weight_46h_Down_BP))
GOData_46h_Down_BP_DE$DEG <- GOData_46h_Down_BP_DE$Significant
GOData_46h_Down_BP_DE$Term <- Term(rownames(GOData_46h_Down_BP_DE))

ggplot(GOData_46h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_46h_Down_BP_DE,
            "topGO_weightFS_46h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 56h Up-regulated

geneList_56h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_56h_Up$Niben))
names(geneList_56h_Up) <- geneNames_Niben261

GOData_56h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_56h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_56h_Up_BP <- runTest(GOData_56h_Up_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_56h_Up_BP <- runTest(GOData_56h_Up_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_56h_Up_BP <- runTest(GOData_56h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_56h_Up_BP <- runTest(GOData_56h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_56h_Up_BP <- usedGO(GOData_56h_Up_BP)
allRes_56h_Up_BP <- GenTable(GOData_56h_Up_BP, 
                             classicFisher = resultsFis_Class_56h_Up_BP,
                             weightFisher = resultsFis_Weight_56h_Up_BP,
                             #classicKS = resultsKS_Class_56h_Up_BP,
                             #weightKS = resultsKS_Weight_56h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_56h_Up_BP))

pVals_Weight_56h_Up_BP <- score(resultsFis_Weight_56h_Up_BP)[score(resultsFis_Weight_56h_Up_BP)<=0.05]

GOData_56h_Up_BP_DE <- termStat(object = GOData_56h_Up_BP, whichGO = names(pVals_Weight_56h_Up_BP))
GOData_56h_Up_BP_DE$DEG <- GOData_56h_Up_BP_DE$Significant
GOData_56h_Up_BP_DE$Term <- Term(rownames(GOData_56h_Up_BP_DE))

ggplot(GOData_56h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_56h_Up_BP_DE,
            "topGO_weightFS_56h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 56h Down-regulated

geneList_56h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_56h_Down$Niben))
names(geneList_56h_Down) <- geneNames_Niben261

GOData_56h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_56h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_56h_Down_BP <- runTest(GOData_56h_Down_BP, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_56h_Down_BP <- runTest(GOData_56h_Down_BP, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_56h_Down_BP <- runTest(GOData_56h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_56h_Down_BP <- runTest(GOData_56h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_56h_Down_BP <- usedGO(GOData_56h_Down_BP)
allRes_56h_Down_BP <- GenTable(GOData_56h_Down_BP, 
                               classicFisher = resultsFis_Class_56h_Down_BP,
                               weightFisher = resultsFis_Weight_56h_Down_BP,
                               #classicKS = resultsKS_Class_56h_Down_BP,
                               #weightKS = resultsKS_Weight_56h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_56h_Down_BP))

pVals_Weight_56h_Down_BP <- score(resultsFis_Weight_56h_Down_BP)[score(resultsFis_Weight_56h_Down_BP)<=0.05]

GOData_56h_Down_BP_DE <- termStat(object = GOData_56h_Down_BP, whichGO = names(pVals_Weight_56h_Down_BP))
GOData_56h_Down_BP_DE$DEG <- GOData_56h_Down_BP_DE$Significant
GOData_56h_Down_BP_DE$Term <- Term(rownames(GOData_56h_Down_BP_DE))

ggplot(GOData_56h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_56h_Down_BP_DE,
            "topGO_weightFS_56h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

################################################################################################################
## UpSetR diagrams of BP in Cornell TimePoints


Cornell_TimePoints_Up_BP <- list("22hUp"= c(GOData_22h_Up_BP_DE$Term),
                                 "25hUp"= c(GOData_25h_Up_BP_DE$Term),
                                 "28hUp"= c(GOData_28h_Up_BP_DE$Term),
                                 "34hUp"= c(GOData_34h_Up_BP_DE$Term),
                                 "37hUp"= c(GOData_37h_Up_BP_DE$Term),
                                 "40hUp"= c(GOData_40h_Up_BP_DE$Term),
                                 "46hUp"= c(GOData_46h_Up_BP_DE$Term),
                                 "56hUp"= c(GOData_56h_Up_BP_DE$Term))

upset(fromList(Cornell_TimePoints_Up_BP), sets= rev(c("22hUp", "25hUp", "28hUp", "34hUp", "37hUp", "40hUp", "46hUp", "56hUp")), order.by = "freq", keep.order = TRUE)

Cornell_TimePoints_Down_BP <- list("22hDown"= c(GOData_22h_Down_BP_DE$Term),
                                 "25hDown"= c(GOData_25h_Down_BP_DE$Term),
                                 "28hDown"= c(GOData_28h_Down_BP_DE$Term),
                                 "34hDown"= c(GOData_34h_Down_BP_DE$Term),
                                 "37hDown"= c(GOData_37h_Down_BP_DE$Term),
                                 "40hDown"= c(GOData_40h_Down_BP_DE$Term),
                                 "46hDown"= c(GOData_46h_Down_BP_DE$Term),
                                 "56hDown"= c(GOData_56h_Down_BP_DE$Term))

upset(fromList(Cornell_TimePoints_Down_BP), sets= rev(c("22hDown", "25hDown", "28hDown", "34hDown", "37hDown", "40hDown", "46hDown", "56hDown")), order.by = "freq", keep.order = TRUE)



#############################  MOLECULAR FUNCTION  ###############################################
#################################################################################################
## 22h Up-regulated

geneList_22h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_22h_Up$Niben))
names(geneList_22h_Up) <- geneNames_Niben261

GOData_22h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_22h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_22h_Up_MF <- runTest(GOData_22h_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_22h_Up_MF <- runTest(GOData_22h_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_22h_Up_MF <- runTest(GOData_22h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_22h_Up_MF <- runTest(GOData_22h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_22h_Up_MF <- usedGO(GOData_22h_Up_MF)
allRes_22h_Up_MF <- GenTable(GOData_22h_Up_MF, 
                             classicFisher = resultsFis_Class_22h_Up_MF,
                             weightFisher = resultsFis_Weight_22h_Up_MF,
                             #classicKS = resultsKS_Class_22h_Up_MF,
                             #weightKS = resultsKS_Weight_22h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_22h_Up_MF))

pVals_Weight_22h_Up_MF <- score(resultsFis_Weight_22h_Up_MF)[score(resultsFis_Weight_22h_Up_MF)<=0.05]

GOData_22h_Up_MF_DE <- termStat(object = GOData_22h_Up_MF, whichGO = names(pVals_Weight_22h_Up_MF))
GOData_22h_Up_MF_DE$DEG <- GOData_22h_Up_MF_DE$Significant
GOData_22h_Up_MF_DE$Term <- Term(rownames(GOData_22h_Up_MF_DE))

ggplot(GOData_22h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_22h_Up_MF_DE,
            "topGO_weightFS_22h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 22h Down-regulated

geneList_22h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_22h_Down$Niben))
names(geneList_22h_Down) <- geneNames_Niben261

GOData_22h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_22h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_22h_Down_MF <- runTest(GOData_22h_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_22h_Down_MF <- runTest(GOData_22h_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_22h_Down_MF <- runTest(GOData_22h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_22h_Down_MF <- runTest(GOData_22h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_22h_Down_MF <- usedGO(GOData_22h_Down_MF)
allRes_22h_Down_MF <- GenTable(GOData_22h_Down_MF, 
                               classicFisher = resultsFis_Class_22h_Down_MF,
                               weightFisher = resultsFis_Weight_22h_Down_MF,
                               #classicKS = resultsKS_Class_22h_Down_MF,
                               #weightKS = resultsKS_Weight_22h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_22h_Down_MF))

pVals_Weight_22h_Down_MF <- score(resultsFis_Weight_22h_Down_MF)[score(resultsFis_Weight_22h_Down_MF)<=0.05]

GOData_22h_Down_MF_DE <- termStat(object = GOData_22h_Down_MF, whichGO = names(pVals_Weight_22h_Down_MF))
GOData_22h_Down_MF_DE$DEG <- GOData_22h_Down_MF_DE$Significant
GOData_22h_Down_MF_DE$Term <- Term(rownames(GOData_22h_Down_MF_DE))

ggplot(GOData_22h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_22h_Down_MF_DE,
            "topGO_weightFS_22h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 25h Up-regulated

geneList_25h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_25h_Up$Niben))
names(geneList_25h_Up) <- geneNames_Niben261

GOData_25h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_25h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_25h_Up_MF <- runTest(GOData_25h_Up_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_25h_Up_MF <- runTest(GOData_25h_Up_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_25h_Up_MF <- runTest(GOData_25h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_25h_Up_MF <- runTest(GOData_25h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_25h_Up_MF <- usedGO(GOData_25h_Up_MF)
allRes_25h_Up_MF <- GenTable(GOData_25h_Up_MF, 
                             classicFisher = resultsFis_Class_25h_Up_MF,
                             weightFisher = resultsFis_Weight_25h_Up_MF,
                             #classicKS = resultsKS_Class_25h_Up_MF,
                             #weightKS = resultsKS_Weight_25h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_25h_Up_MF))

pVals_Weight_25h_Up_MF <- score(resultsFis_Weight_25h_Up_MF)[score(resultsFis_Weight_25h_Up_MF)<=0.05]

GOData_25h_Up_MF_DE <- termStat(object = GOData_25h_Up_MF, whichGO = names(pVals_Weight_25h_Up_MF))
GOData_25h_Up_MF_DE$DEG <- GOData_25h_Up_MF_DE$Significant
GOData_25h_Up_MF_DE$Term <- Term(rownames(GOData_25h_Up_MF_DE))

ggplot(GOData_25h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_25h_Up_MF_DE,
            "topGO_weightFS_25h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 25h Down-regulated

geneList_25h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_25h_Down$Niben))
names(geneList_25h_Down) <- geneNames_Niben261

GOData_25h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_25h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_25h_Down_MF <- runTest(GOData_25h_Down_MF, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_25h_Down_MF <- runTest(GOData_25h_Down_MF, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_25h_Down_MF <- runTest(GOData_25h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_25h_Down_MF <- runTest(GOData_25h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_25h_Down_MF <- usedGO(GOData_25h_Down_MF)
allRes_25h_Down_MF <- GenTable(GOData_25h_Down_MF, 
                               classicFisher = resultsFis_Class_25h_Down_MF,
                               weightFisher = resultsFis_Weight_25h_Down_MF,
                               #classicKS = resultsKS_Class_25h_Down_MF,
                               #weightKS = resultsKS_Weight_25h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_25h_Down_MF))

pVals_Weight_25h_Down_MF <- score(resultsFis_Weight_25h_Down_MF)[score(resultsFis_Weight_25h_Down_MF)<=0.05]

GOData_25h_Down_MF_DE <- termStat(object = GOData_25h_Down_MF, whichGO = names(pVals_Weight_25h_Down_MF))
GOData_25h_Down_MF_DE$DEG <- GOData_25h_Down_MF_DE$Significant
GOData_25h_Down_MF_DE$Term <- Term(rownames(GOData_25h_Down_MF_DE))

ggplot(GOData_25h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_25h_Down_MF_DE,
            "topGO_weightFS_25h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 28h Up-regulated

geneList_28h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_28h_Up$Niben))
names(geneList_28h_Up) <- geneNames_Niben261

GOData_28h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_28h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_28h_Up_MF <- runTest(GOData_28h_Up_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_28h_Up_MF <- runTest(GOData_28h_Up_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_28h_Up_MF <- runTest(GOData_28h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_28h_Up_MF <- runTest(GOData_28h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_28h_Up_MF <- usedGO(GOData_28h_Up_MF)
allRes_28h_Up_MF <- GenTable(GOData_28h_Up_MF, 
                             classicFisher = resultsFis_Class_28h_Up_MF,
                             weightFisher = resultsFis_Weight_28h_Up_MF,
                             #classicKS = resultsKS_Class_28h_Up_MF,
                             #weightKS = resultsKS_Weight_28h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_28h_Up_MF))

pVals_Weight_28h_Up_MF <- score(resultsFis_Weight_28h_Up_MF)[score(resultsFis_Weight_28h_Up_MF)<=0.05]

GOData_28h_Up_MF_DE <- termStat(object = GOData_28h_Up_MF, whichGO = names(pVals_Weight_28h_Up_MF))
GOData_28h_Up_MF_DE$DEG <- GOData_28h_Up_MF_DE$Significant
GOData_28h_Up_MF_DE$Term <- Term(rownames(GOData_28h_Up_MF_DE))

ggplot(GOData_28h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_28h_Up_MF_DE,
            "topGO_weightFS_28h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 28h Down-regulated

geneList_28h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_28h_Down$Niben))
names(geneList_28h_Down) <- geneNames_Niben261

GOData_28h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_28h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_28h_Down_MF <- runTest(GOData_28h_Down_MF, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_28h_Down_MF <- runTest(GOData_28h_Down_MF, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_28h_Down_MF <- runTest(GOData_28h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_28h_Down_MF <- runTest(GOData_28h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_28h_Down_MF <- usedGO(GOData_28h_Down_MF)
allRes_28h_Down_MF <- GenTable(GOData_28h_Down_MF, 
                               classicFisher = resultsFis_Class_28h_Down_MF,
                               weightFisher = resultsFis_Weight_28h_Down_MF,
                               #classicKS = resultsKS_Class_28h_Down_MF,
                               #weightKS = resultsKS_Weight_28h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_28h_Down_MF))

pVals_Weight_28h_Down_MF <- score(resultsFis_Weight_28h_Down_MF)[score(resultsFis_Weight_28h_Down_MF)<=0.05]

GOData_28h_Down_MF_DE <- termStat(object = GOData_28h_Down_MF, whichGO = names(pVals_Weight_28h_Down_MF))
GOData_28h_Down_MF_DE$DEG <- GOData_28h_Down_MF_DE$Significant
GOData_28h_Down_MF_DE$Term <- Term(rownames(GOData_28h_Down_MF_DE))

ggplot(GOData_28h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_28h_Down_MF_DE,
            "topGO_weightFS_28h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 34h Up-regulated

geneList_34h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_34h_Up$Niben))
names(geneList_34h_Up) <- geneNames_Niben261

GOData_34h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_34h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_34h_Up_MF <- runTest(GOData_34h_Up_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_34h_Up_MF <- runTest(GOData_34h_Up_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_34h_Up_MF <- runTest(GOData_34h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_34h_Up_MF <- runTest(GOData_34h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_34h_Up_MF <- usedGO(GOData_34h_Up_MF)
allRes_34h_Up_MF <- GenTable(GOData_34h_Up_MF, 
                             classicFisher = resultsFis_Class_34h_Up_MF,
                             weightFisher = resultsFis_Weight_34h_Up_MF,
                             #classicKS = resultsKS_Class_34h_Up_MF,
                             #weightKS = resultsKS_Weight_34h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_34h_Up_MF))

pVals_Weight_34h_Up_MF <- score(resultsFis_Weight_34h_Up_MF)[score(resultsFis_Weight_34h_Up_MF)<=0.05]

GOData_34h_Up_MF_DE <- termStat(object = GOData_34h_Up_MF, whichGO = names(pVals_Weight_34h_Up_MF))
GOData_34h_Up_MF_DE$DEG <- GOData_34h_Up_MF_DE$Significant
GOData_34h_Up_MF_DE$Term <- Term(rownames(GOData_34h_Up_MF_DE))

ggplot(GOData_34h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_34h_Up_MF_DE,
            "topGO_weightFS_34h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 34h Down-regulated

geneList_34h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_34h_Down$Niben))
names(geneList_34h_Down) <- geneNames_Niben261

GOData_34h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_34h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_34h_Down_MF <- runTest(GOData_34h_Down_MF, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_34h_Down_MF <- runTest(GOData_34h_Down_MF, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_34h_Down_MF <- runTest(GOData_34h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_34h_Down_MF <- runTest(GOData_34h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_34h_Down_MF <- usedGO(GOData_34h_Down_MF)
allRes_34h_Down_MF <- GenTable(GOData_34h_Down_MF, 
                               classicFisher = resultsFis_Class_34h_Down_MF,
                               weightFisher = resultsFis_Weight_34h_Down_MF,
                               #classicKS = resultsKS_Class_34h_Down_MF,
                               #weightKS = resultsKS_Weight_34h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_34h_Down_MF))

pVals_Weight_34h_Down_MF <- score(resultsFis_Weight_34h_Down_MF)[score(resultsFis_Weight_34h_Down_MF)<=0.05]

GOData_34h_Down_MF_DE <- termStat(object = GOData_34h_Down_MF, whichGO = names(pVals_Weight_34h_Down_MF))
GOData_34h_Down_MF_DE$DEG <- GOData_34h_Down_MF_DE$Significant
GOData_34h_Down_MF_DE$Term <- Term(rownames(GOData_34h_Down_MF_DE))

ggplot(GOData_34h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_34h_Down_MF_DE,
            "topGO_weightFS_34h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 37h Up-regulated

geneList_37h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_37h_Up$Niben))
names(geneList_37h_Up) <- geneNames_Niben261

GOData_37h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_37h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_37h_Up_MF <- runTest(GOData_37h_Up_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_37h_Up_MF <- runTest(GOData_37h_Up_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_37h_Up_MF <- runTest(GOData_37h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_37h_Up_MF <- runTest(GOData_37h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_37h_Up_MF <- usedGO(GOData_37h_Up_MF)
allRes_37h_Up_MF <- GenTable(GOData_37h_Up_MF, 
                             classicFisher = resultsFis_Class_37h_Up_MF,
                             weightFisher = resultsFis_Weight_37h_Up_MF,
                             #classicKS = resultsKS_Class_37h_Up_MF,
                             #weightKS = resultsKS_Weight_37h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_37h_Up_MF))

pVals_Weight_37h_Up_MF <- score(resultsFis_Weight_37h_Up_MF)[score(resultsFis_Weight_37h_Up_MF)<=0.05]

GOData_37h_Up_MF_DE <- termStat(object = GOData_37h_Up_MF, whichGO = names(pVals_Weight_37h_Up_MF))
GOData_37h_Up_MF_DE$DEG <- GOData_37h_Up_MF_DE$Significant
GOData_37h_Up_MF_DE$Term <- Term(rownames(GOData_37h_Up_MF_DE))

ggplot(GOData_37h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_37h_Up_MF_DE,
            "topGO_weightFS_37h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 37h Down-regulated

geneList_37h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_37h_Down$Niben))
names(geneList_37h_Down) <- geneNames_Niben261

GOData_37h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_37h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_37h_Down_MF <- runTest(GOData_37h_Down_MF, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_37h_Down_MF <- runTest(GOData_37h_Down_MF, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_37h_Down_MF <- runTest(GOData_37h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_37h_Down_MF <- runTest(GOData_37h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_37h_Down_MF <- usedGO(GOData_37h_Down_MF)
allRes_37h_Down_MF <- GenTable(GOData_37h_Down_MF, 
                               classicFisher = resultsFis_Class_37h_Down_MF,
                               weightFisher = resultsFis_Weight_37h_Down_MF,
                               #classicKS = resultsKS_Class_37h_Down_MF,
                               #weightKS = resultsKS_Weight_37h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_37h_Down_MF))

pVals_Weight_37h_Down_MF <- score(resultsFis_Weight_37h_Down_MF)[score(resultsFis_Weight_37h_Down_MF)<=0.05]

GOData_37h_Down_MF_DE <- termStat(object = GOData_37h_Down_MF, whichGO = names(pVals_Weight_37h_Down_MF))
GOData_37h_Down_MF_DE$DEG <- GOData_37h_Down_MF_DE$Significant
GOData_37h_Down_MF_DE$Term <- Term(rownames(GOData_37h_Down_MF_DE))

ggplot(GOData_37h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_37h_Down_MF_DE,
            "topGO_weightFS_37h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 40h Up-regulated

geneList_40h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_40h_Up$Niben))
names(geneList_40h_Up) <- geneNames_Niben261

GOData_40h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_40h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_40h_Up_MF <- runTest(GOData_40h_Up_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_40h_Up_MF <- runTest(GOData_40h_Up_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_40h_Up_MF <- runTest(GOData_40h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_40h_Up_MF <- runTest(GOData_40h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_40h_Up_MF <- usedGO(GOData_40h_Up_MF)
allRes_40h_Up_MF <- GenTable(GOData_40h_Up_MF, 
                             classicFisher = resultsFis_Class_40h_Up_MF,
                             weightFisher = resultsFis_Weight_40h_Up_MF,
                             #classicKS = resultsKS_Class_40h_Up_MF,
                             #weightKS = resultsKS_Weight_40h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_40h_Up_MF))

pVals_Weight_40h_Up_MF <- score(resultsFis_Weight_40h_Up_MF)[score(resultsFis_Weight_40h_Up_MF)<=0.05]

GOData_40h_Up_MF_DE <- termStat(object = GOData_40h_Up_MF, whichGO = names(pVals_Weight_40h_Up_MF))
GOData_40h_Up_MF_DE$DEG <- GOData_40h_Up_MF_DE$Significant
GOData_40h_Up_MF_DE$Term <- Term(rownames(GOData_40h_Up_MF_DE))

ggplot(GOData_40h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_40h_Up_MF_DE,
            "topGO_weightFS_40h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 40h Down-regulated

geneList_40h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_40h_Down$Niben))
names(geneList_40h_Down) <- geneNames_Niben261

GOData_40h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_40h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_40h_Down_MF <- runTest(GOData_40h_Down_MF, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_40h_Down_MF <- runTest(GOData_40h_Down_MF, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_40h_Down_MF <- runTest(GOData_40h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_40h_Down_MF <- runTest(GOData_40h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_40h_Down_MF <- usedGO(GOData_40h_Down_MF)
allRes_40h_Down_MF <- GenTable(GOData_40h_Down_MF, 
                               classicFisher = resultsFis_Class_40h_Down_MF,
                               weightFisher = resultsFis_Weight_40h_Down_MF,
                               #classicKS = resultsKS_Class_40h_Down_MF,
                               #weightKS = resultsKS_Weight_40h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_40h_Down_MF))

pVals_Weight_40h_Down_MF <- score(resultsFis_Weight_40h_Down_MF)[score(resultsFis_Weight_40h_Down_MF)<=0.05]

GOData_40h_Down_MF_DE <- termStat(object = GOData_40h_Down_MF, whichGO = names(pVals_Weight_40h_Down_MF))
GOData_40h_Down_MF_DE$DEG <- GOData_40h_Down_MF_DE$Significant
GOData_40h_Down_MF_DE$Term <- Term(rownames(GOData_40h_Down_MF_DE))

ggplot(GOData_40h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_40h_Down_MF_DE,
            "topGO_weightFS_40h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 46h Up-regulated

geneList_46h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_46h_Up$Niben))
names(geneList_46h_Up) <- geneNames_Niben261

GOData_46h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_46h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_46h_Up_MF <- runTest(GOData_46h_Up_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_46h_Up_MF <- runTest(GOData_46h_Up_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_46h_Up_MF <- runTest(GOData_46h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_46h_Up_MF <- runTest(GOData_46h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_46h_Up_MF <- usedGO(GOData_46h_Up_MF)
allRes_46h_Up_MF <- GenTable(GOData_46h_Up_MF, 
                             classicFisher = resultsFis_Class_46h_Up_MF,
                             weightFisher = resultsFis_Weight_46h_Up_MF,
                             #classicKS = resultsKS_Class_46h_Up_MF,
                             #weightKS = resultsKS_Weight_46h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_46h_Up_MF))

pVals_Weight_46h_Up_MF <- score(resultsFis_Weight_46h_Up_MF)[score(resultsFis_Weight_46h_Up_MF)<=0.05]

GOData_46h_Up_MF_DE <- termStat(object = GOData_46h_Up_MF, whichGO = names(pVals_Weight_46h_Up_MF))
GOData_46h_Up_MF_DE$DEG <- GOData_46h_Up_MF_DE$Significant
GOData_46h_Up_MF_DE$Term <- Term(rownames(GOData_46h_Up_MF_DE))

ggplot(GOData_46h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_46h_Up_MF_DE,
            "topGO_weightFS_46h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 46h Down-regulated

geneList_46h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_46h_Down$Niben))
names(geneList_46h_Down) <- geneNames_Niben261

GOData_46h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_46h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_46h_Down_MF <- runTest(GOData_46h_Down_MF, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_46h_Down_MF <- runTest(GOData_46h_Down_MF, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_46h_Down_MF <- runTest(GOData_46h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_46h_Down_MF <- runTest(GOData_46h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_46h_Down_MF <- usedGO(GOData_46h_Down_MF)
allRes_46h_Down_MF <- GenTable(GOData_46h_Down_MF, 
                               classicFisher = resultsFis_Class_46h_Down_MF,
                               weightFisher = resultsFis_Weight_46h_Down_MF,
                               #classicKS = resultsKS_Class_46h_Down_MF,
                               #weightKS = resultsKS_Weight_46h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_46h_Down_MF))

pVals_Weight_46h_Down_MF <- score(resultsFis_Weight_46h_Down_MF)[score(resultsFis_Weight_46h_Down_MF)<=0.05]

GOData_46h_Down_MF_DE <- termStat(object = GOData_46h_Down_MF, whichGO = names(pVals_Weight_46h_Down_MF))
GOData_46h_Down_MF_DE$DEG <- GOData_46h_Down_MF_DE$Significant
GOData_46h_Down_MF_DE$Term <- Term(rownames(GOData_46h_Down_MF_DE))

ggplot(GOData_46h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_46h_Down_MF_DE,
            "topGO_weightFS_46h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 56h Up-regulated

geneList_56h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_56h_Up$Niben))
names(geneList_56h_Up) <- geneNames_Niben261

GOData_56h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_56h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_56h_Up_MF <- runTest(GOData_56h_Up_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_56h_Up_MF <- runTest(GOData_56h_Up_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_56h_Up_MF <- runTest(GOData_56h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_56h_Up_MF <- runTest(GOData_56h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_56h_Up_MF <- usedGO(GOData_56h_Up_MF)
allRes_56h_Up_MF <- GenTable(GOData_56h_Up_MF, 
                             classicFisher = resultsFis_Class_56h_Up_MF,
                             weightFisher = resultsFis_Weight_56h_Up_MF,
                             #classicKS = resultsKS_Class_56h_Up_MF,
                             #weightKS = resultsKS_Weight_56h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_56h_Up_MF))

pVals_Weight_56h_Up_MF <- score(resultsFis_Weight_56h_Up_MF)[score(resultsFis_Weight_56h_Up_MF)<=0.05]

GOData_56h_Up_MF_DE <- termStat(object = GOData_56h_Up_MF, whichGO = names(pVals_Weight_56h_Up_MF))
GOData_56h_Up_MF_DE$DEG <- GOData_56h_Up_MF_DE$Significant
GOData_56h_Up_MF_DE$Term <- Term(rownames(GOData_56h_Up_MF_DE))

ggplot(GOData_56h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_56h_Up_MF_DE,
            "topGO_weightFS_56h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 56h Down-regulated

geneList_56h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_56h_Down$Niben))
names(geneList_56h_Down) <- geneNames_Niben261

GOData_56h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_56h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_56h_Down_MF <- runTest(GOData_56h_Down_MF, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_56h_Down_MF <- runTest(GOData_56h_Down_MF, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_56h_Down_MF <- runTest(GOData_56h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_56h_Down_MF <- runTest(GOData_56h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_56h_Down_MF <- usedGO(GOData_56h_Down_MF)
allRes_56h_Down_MF <- GenTable(GOData_56h_Down_MF, 
                               classicFisher = resultsFis_Class_56h_Down_MF,
                               weightFisher = resultsFis_Weight_56h_Down_MF,
                               #classicKS = resultsKS_Class_56h_Down_MF,
                               #weightKS = resultsKS_Weight_56h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_56h_Down_MF))

pVals_Weight_56h_Down_MF <- score(resultsFis_Weight_56h_Down_MF)[score(resultsFis_Weight_56h_Down_MF)<=0.05]

GOData_56h_Down_MF_DE <- termStat(object = GOData_56h_Down_MF, whichGO = names(pVals_Weight_56h_Down_MF))
GOData_56h_Down_MF_DE$DEG <- GOData_56h_Down_MF_DE$Significant
GOData_56h_Down_MF_DE$Term <- Term(rownames(GOData_56h_Down_MF_DE))

ggplot(GOData_56h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_56h_Down_MF_DE,
            "topGO_weightFS_56h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

################################################################################################################
## UpSetR diagrams of MF in Cornell TimePoints


Cornell_TimePoints_Up_MF <- list("22hUp"= c(GOData_22h_Up_MF_DE$Term),
                                 "25hUp"= c(GOData_25h_Up_MF_DE$Term),
                                 "28hUp"= c(GOData_28h_Up_MF_DE$Term),
                                 "34hUp"= c(GOData_34h_Up_MF_DE$Term),
                                 "37hUp"= c(GOData_37h_Up_MF_DE$Term),
                                 "40hUp"= c(GOData_40h_Up_MF_DE$Term),
                                 "46hUp"= c(GOData_46h_Up_MF_DE$Term),
                                 "56hUp"= c(GOData_56h_Up_MF_DE$Term))

upset(fromList(Cornell_TimePoints_Up_MF), sets= rev(c("22hUp", "25hUp", "28hUp", "34hUp", "37hUp", "40hUp", "46hUp", "56hUp")), order.by = "freq", keep.order = TRUE)

Cornell_TimePoints_Down_MF <- list("22hDown"= c(GOData_22h_Down_MF_DE$Term),
                                   "25hDown"= c(GOData_25h_Down_MF_DE$Term),
                                   "28hDown"= c(GOData_28h_Down_MF_DE$Term),
                                   "34hDown"= c(GOData_34h_Down_MF_DE$Term),
                                   "37hDown"= c(GOData_37h_Down_MF_DE$Term),
                                   "40hDown"= c(GOData_40h_Down_MF_DE$Term),
                                   "46hDown"= c(GOData_46h_Down_MF_DE$Term),
                                   "56hDown"= c(GOData_56h_Down_MF_DE$Term))

upset(fromList(Cornell_TimePoints_Down_MF), sets= rev(c("22hDown", "25hDown", "28hDown", "34hDown", "37hDown", "40hDown", "46hDown", "56hDown")), order.by = "freq", keep.order = TRUE)



#############################  CELLULAR COMPONENT  ###############################################
#################################################################################################
## 22h Up-regulated

geneList_22h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_22h_Up$Niben))
names(geneList_22h_Up) <- geneNames_Niben261

GOData_22h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_22h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_22h_Up_CC <- runTest(GOData_22h_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_22h_Up_CC <- runTest(GOData_22h_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_22h_Up_CC <- runTest(GOData_22h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_22h_Up_CC <- runTest(GOData_22h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_22h_Up_CC <- usedGO(GOData_22h_Up_CC)
allRes_22h_Up_CC <- GenTable(GOData_22h_Up_CC, 
                             classicFisher = resultsFis_Class_22h_Up_CC,
                             weightFisher = resultsFis_Weight_22h_Up_CC,
                             #classicKS = resultsKS_Class_22h_Up_CC,
                             #weightKS = resultsKS_Weight_22h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_22h_Up_CC))

pVals_Weight_22h_Up_CC <- score(resultsFis_Weight_22h_Up_CC)[score(resultsFis_Weight_22h_Up_CC)<=0.05]

GOData_22h_Up_CC_DE <- termStat(object = GOData_22h_Up_CC, whichGO = names(pVals_Weight_22h_Up_CC))
GOData_22h_Up_CC_DE$DEG <- GOData_22h_Up_CC_DE$Significant
GOData_22h_Up_CC_DE$Term <- Term(rownames(GOData_22h_Up_CC_DE))

ggplot(GOData_22h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_22h_Up_CC_DE,
            "topGO_weightFS_22h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 22h Down-regulated

geneList_22h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_22h_Down$Niben))
names(geneList_22h_Down) <- geneNames_Niben261

GOData_22h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_22h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_22h_Down_CC <- runTest(GOData_22h_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_22h_Down_CC <- runTest(GOData_22h_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_22h_Down_CC <- runTest(GOData_22h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_22h_Down_CC <- runTest(GOData_22h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_22h_Down_CC <- usedGO(GOData_22h_Down_CC)
allRes_22h_Down_CC <- GenTable(GOData_22h_Down_CC, 
                               classicFisher = resultsFis_Class_22h_Down_CC,
                               weightFisher = resultsFis_Weight_22h_Down_CC,
                               #classicKS = resultsKS_Class_22h_Down_CC,
                               #weightKS = resultsKS_Weight_22h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_22h_Down_CC))

pVals_Weight_22h_Down_CC <- score(resultsFis_Weight_22h_Down_CC)[score(resultsFis_Weight_22h_Down_CC)<=0.05]

GOData_22h_Down_CC_DE <- termStat(object = GOData_22h_Down_CC, whichGO = names(pVals_Weight_22h_Down_CC))
GOData_22h_Down_CC_DE$DEG <- GOData_22h_Down_CC_DE$Significant
GOData_22h_Down_CC_DE$Term <- Term(rownames(GOData_22h_Down_CC_DE))

ggplot(GOData_22h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_22h_Down_CC_DE,
            "topGO_weightFS_22h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 25h Up-regulated

geneList_25h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_25h_Up$Niben))
names(geneList_25h_Up) <- geneNames_Niben261

GOData_25h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_25h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_25h_Up_CC <- runTest(GOData_25h_Up_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_25h_Up_CC <- runTest(GOData_25h_Up_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_25h_Up_CC <- runTest(GOData_25h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_25h_Up_CC <- runTest(GOData_25h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_25h_Up_CC <- usedGO(GOData_25h_Up_CC)
allRes_25h_Up_CC <- GenTable(GOData_25h_Up_CC, 
                             classicFisher = resultsFis_Class_25h_Up_CC,
                             weightFisher = resultsFis_Weight_25h_Up_CC,
                             #classicKS = resultsKS_Class_25h_Up_CC,
                             #weightKS = resultsKS_Weight_25h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_25h_Up_CC))

pVals_Weight_25h_Up_CC <- score(resultsFis_Weight_25h_Up_CC)[score(resultsFis_Weight_25h_Up_CC)<=0.05]

GOData_25h_Up_CC_DE <- termStat(object = GOData_25h_Up_CC, whichGO = names(pVals_Weight_25h_Up_CC))
GOData_25h_Up_CC_DE$DEG <- GOData_25h_Up_CC_DE$Significant
GOData_25h_Up_CC_DE$Term <- Term(rownames(GOData_25h_Up_CC_DE))

ggplot(GOData_25h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_25h_Up_CC_DE,
            "topGO_weightFS_25h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 25h Down-regulated

geneList_25h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_25h_Down$Niben))
names(geneList_25h_Down) <- geneNames_Niben261

GOData_25h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_25h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_25h_Down_CC <- runTest(GOData_25h_Down_CC, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_25h_Down_CC <- runTest(GOData_25h_Down_CC, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_25h_Down_CC <- runTest(GOData_25h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_25h_Down_CC <- runTest(GOData_25h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_25h_Down_CC <- usedGO(GOData_25h_Down_CC)
allRes_25h_Down_CC <- GenTable(GOData_25h_Down_CC, 
                               classicFisher = resultsFis_Class_25h_Down_CC,
                               weightFisher = resultsFis_Weight_25h_Down_CC,
                               #classicKS = resultsKS_Class_25h_Down_CC,
                               #weightKS = resultsKS_Weight_25h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_25h_Down_CC))

pVals_Weight_25h_Down_CC <- score(resultsFis_Weight_25h_Down_CC)[score(resultsFis_Weight_25h_Down_CC)<=0.05]

GOData_25h_Down_CC_DE <- termStat(object = GOData_25h_Down_CC, whichGO = names(pVals_Weight_25h_Down_CC))
GOData_25h_Down_CC_DE$DEG <- GOData_25h_Down_CC_DE$Significant
GOData_25h_Down_CC_DE$Term <- Term(rownames(GOData_25h_Down_CC_DE))

ggplot(GOData_25h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_25h_Down_CC_DE,
            "topGO_weightFS_25h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 28h Up-regulated

geneList_28h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_28h_Up$Niben))
names(geneList_28h_Up) <- geneNames_Niben261

GOData_28h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_28h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_28h_Up_CC <- runTest(GOData_28h_Up_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_28h_Up_CC <- runTest(GOData_28h_Up_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_28h_Up_CC <- runTest(GOData_28h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_28h_Up_CC <- runTest(GOData_28h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_28h_Up_CC <- usedGO(GOData_28h_Up_CC)
allRes_28h_Up_CC <- GenTable(GOData_28h_Up_CC, 
                             classicFisher = resultsFis_Class_28h_Up_CC,
                             weightFisher = resultsFis_Weight_28h_Up_CC,
                             #classicKS = resultsKS_Class_28h_Up_CC,
                             #weightKS = resultsKS_Weight_28h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_28h_Up_CC))

pVals_Weight_28h_Up_CC <- score(resultsFis_Weight_28h_Up_CC)[score(resultsFis_Weight_28h_Up_CC)<=0.05]

GOData_28h_Up_CC_DE <- termStat(object = GOData_28h_Up_CC, whichGO = names(pVals_Weight_28h_Up_CC))
GOData_28h_Up_CC_DE$DEG <- GOData_28h_Up_CC_DE$Significant
GOData_28h_Up_CC_DE$Term <- Term(rownames(GOData_28h_Up_CC_DE))

ggplot(GOData_28h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_28h_Up_CC_DE,
            "topGO_weightFS_28h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 28h Down-regulated

geneList_28h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_28h_Down$Niben))
names(geneList_28h_Down) <- geneNames_Niben261

GOData_28h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_28h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_28h_Down_CC <- runTest(GOData_28h_Down_CC, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_28h_Down_CC <- runTest(GOData_28h_Down_CC, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_28h_Down_CC <- runTest(GOData_28h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_28h_Down_CC <- runTest(GOData_28h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_28h_Down_CC <- usedGO(GOData_28h_Down_CC)
allRes_28h_Down_CC <- GenTable(GOData_28h_Down_CC, 
                               classicFisher = resultsFis_Class_28h_Down_CC,
                               weightFisher = resultsFis_Weight_28h_Down_CC,
                               #classicKS = resultsKS_Class_28h_Down_CC,
                               #weightKS = resultsKS_Weight_28h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_28h_Down_CC))

pVals_Weight_28h_Down_CC <- score(resultsFis_Weight_28h_Down_CC)[score(resultsFis_Weight_28h_Down_CC)<=0.05]

GOData_28h_Down_CC_DE <- termStat(object = GOData_28h_Down_CC, whichGO = names(pVals_Weight_28h_Down_CC))
GOData_28h_Down_CC_DE$DEG <- GOData_28h_Down_CC_DE$Significant
GOData_28h_Down_CC_DE$Term <- Term(rownames(GOData_28h_Down_CC_DE))

ggplot(GOData_28h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_28h_Down_CC_DE,
            "topGO_weightFS_28h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 34h Up-regulated

geneList_34h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_34h_Up$Niben))
names(geneList_34h_Up) <- geneNames_Niben261

GOData_34h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_34h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_34h_Up_CC <- runTest(GOData_34h_Up_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_34h_Up_CC <- runTest(GOData_34h_Up_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_34h_Up_CC <- runTest(GOData_34h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_34h_Up_CC <- runTest(GOData_34h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_34h_Up_CC <- usedGO(GOData_34h_Up_CC)
allRes_34h_Up_CC <- GenTable(GOData_34h_Up_CC, 
                             classicFisher = resultsFis_Class_34h_Up_CC,
                             weightFisher = resultsFis_Weight_34h_Up_CC,
                             #classicKS = resultsKS_Class_34h_Up_CC,
                             #weightKS = resultsKS_Weight_34h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_34h_Up_CC))

pVals_Weight_34h_Up_CC <- score(resultsFis_Weight_34h_Up_CC)[score(resultsFis_Weight_34h_Up_CC)<=0.05]

GOData_34h_Up_CC_DE <- termStat(object = GOData_34h_Up_CC, whichGO = names(pVals_Weight_34h_Up_CC))
GOData_34h_Up_CC_DE$DEG <- GOData_34h_Up_CC_DE$Significant
GOData_34h_Up_CC_DE$Term <- Term(rownames(GOData_34h_Up_CC_DE))

ggplot(GOData_34h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_34h_Up_CC_DE,
            "topGO_weightFS_34h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 34h Down-regulated

geneList_34h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_34h_Down$Niben))
names(geneList_34h_Down) <- geneNames_Niben261

GOData_34h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_34h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_34h_Down_CC <- runTest(GOData_34h_Down_CC, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_34h_Down_CC <- runTest(GOData_34h_Down_CC, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_34h_Down_CC <- runTest(GOData_34h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_34h_Down_CC <- runTest(GOData_34h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_34h_Down_CC <- usedGO(GOData_34h_Down_CC)
allRes_34h_Down_CC <- GenTable(GOData_34h_Down_CC, 
                               classicFisher = resultsFis_Class_34h_Down_CC,
                               weightFisher = resultsFis_Weight_34h_Down_CC,
                               #classicKS = resultsKS_Class_34h_Down_CC,
                               #weightKS = resultsKS_Weight_34h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_34h_Down_CC))

pVals_Weight_34h_Down_CC <- score(resultsFis_Weight_34h_Down_CC)[score(resultsFis_Weight_34h_Down_CC)<=0.05]

GOData_34h_Down_CC_DE <- termStat(object = GOData_34h_Down_CC, whichGO = names(pVals_Weight_34h_Down_CC))
GOData_34h_Down_CC_DE$DEG <- GOData_34h_Down_CC_DE$Significant
GOData_34h_Down_CC_DE$Term <- Term(rownames(GOData_34h_Down_CC_DE))

ggplot(GOData_34h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_34h_Down_CC_DE,
            "topGO_weightFS_34h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 37h Up-regulated

geneList_37h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_37h_Up$Niben))
names(geneList_37h_Up) <- geneNames_Niben261

GOData_37h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_37h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_37h_Up_CC <- runTest(GOData_37h_Up_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_37h_Up_CC <- runTest(GOData_37h_Up_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_37h_Up_CC <- runTest(GOData_37h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_37h_Up_CC <- runTest(GOData_37h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_37h_Up_CC <- usedGO(GOData_37h_Up_CC)
allRes_37h_Up_CC <- GenTable(GOData_37h_Up_CC, 
                             classicFisher = resultsFis_Class_37h_Up_CC,
                             weightFisher = resultsFis_Weight_37h_Up_CC,
                             #classicKS = resultsKS_Class_37h_Up_CC,
                             #weightKS = resultsKS_Weight_37h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_37h_Up_CC))

pVals_Weight_37h_Up_CC <- score(resultsFis_Weight_37h_Up_CC)[score(resultsFis_Weight_37h_Up_CC)<=0.05]

GOData_37h_Up_CC_DE <- termStat(object = GOData_37h_Up_CC, whichGO = names(pVals_Weight_37h_Up_CC))
GOData_37h_Up_CC_DE$DEG <- GOData_37h_Up_CC_DE$Significant
GOData_37h_Up_CC_DE$Term <- Term(rownames(GOData_37h_Up_CC_DE))

ggplot(GOData_37h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_37h_Up_CC_DE,
            "topGO_weightFS_37h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 37h Down-regulated

geneList_37h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_37h_Down$Niben))
names(geneList_37h_Down) <- geneNames_Niben261

GOData_37h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_37h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_37h_Down_CC <- runTest(GOData_37h_Down_CC, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_37h_Down_CC <- runTest(GOData_37h_Down_CC, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_37h_Down_CC <- runTest(GOData_37h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_37h_Down_CC <- runTest(GOData_37h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_37h_Down_CC <- usedGO(GOData_37h_Down_CC)
allRes_37h_Down_CC <- GenTable(GOData_37h_Down_CC, 
                               classicFisher = resultsFis_Class_37h_Down_CC,
                               weightFisher = resultsFis_Weight_37h_Down_CC,
                               #classicKS = resultsKS_Class_37h_Down_CC,
                               #weightKS = resultsKS_Weight_37h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_37h_Down_CC))

pVals_Weight_37h_Down_CC <- score(resultsFis_Weight_37h_Down_CC)[score(resultsFis_Weight_37h_Down_CC)<=0.05]

GOData_37h_Down_CC_DE <- termStat(object = GOData_37h_Down_CC, whichGO = names(pVals_Weight_37h_Down_CC))
GOData_37h_Down_CC_DE$DEG <- GOData_37h_Down_CC_DE$Significant
GOData_37h_Down_CC_DE$Term <- Term(rownames(GOData_37h_Down_CC_DE))

ggplot(GOData_37h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_37h_Down_CC_DE,
            "topGO_weightFS_37h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 40h Up-regulated

geneList_40h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_40h_Up$Niben))
names(geneList_40h_Up) <- geneNames_Niben261

GOData_40h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_40h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_40h_Up_CC <- runTest(GOData_40h_Up_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_40h_Up_CC <- runTest(GOData_40h_Up_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_40h_Up_CC <- runTest(GOData_40h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_40h_Up_CC <- runTest(GOData_40h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_40h_Up_CC <- usedGO(GOData_40h_Up_CC)
allRes_40h_Up_CC <- GenTable(GOData_40h_Up_CC, 
                             classicFisher = resultsFis_Class_40h_Up_CC,
                             weightFisher = resultsFis_Weight_40h_Up_CC,
                             #classicKS = resultsKS_Class_40h_Up_CC,
                             #weightKS = resultsKS_Weight_40h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_40h_Up_CC))

pVals_Weight_40h_Up_CC <- score(resultsFis_Weight_40h_Up_CC)[score(resultsFis_Weight_40h_Up_CC)<=0.05]

GOData_40h_Up_CC_DE <- termStat(object = GOData_40h_Up_CC, whichGO = names(pVals_Weight_40h_Up_CC))
GOData_40h_Up_CC_DE$DEG <- GOData_40h_Up_CC_DE$Significant
GOData_40h_Up_CC_DE$Term <- Term(rownames(GOData_40h_Up_CC_DE))

ggplot(GOData_40h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_40h_Up_CC_DE,
            "topGO_weightFS_40h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 40h Down-regulated

geneList_40h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_40h_Down$Niben))
names(geneList_40h_Down) <- geneNames_Niben261

GOData_40h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_40h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_40h_Down_CC <- runTest(GOData_40h_Down_CC, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_40h_Down_CC <- runTest(GOData_40h_Down_CC, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_40h_Down_CC <- runTest(GOData_40h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_40h_Down_CC <- runTest(GOData_40h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_40h_Down_CC <- usedGO(GOData_40h_Down_CC)
allRes_40h_Down_CC <- GenTable(GOData_40h_Down_CC, 
                               classicFisher = resultsFis_Class_40h_Down_CC,
                               weightFisher = resultsFis_Weight_40h_Down_CC,
                               #classicKS = resultsKS_Class_40h_Down_CC,
                               #weightKS = resultsKS_Weight_40h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_40h_Down_CC))

pVals_Weight_40h_Down_CC <- score(resultsFis_Weight_40h_Down_CC)[score(resultsFis_Weight_40h_Down_CC)<=0.05]

GOData_40h_Down_CC_DE <- termStat(object = GOData_40h_Down_CC, whichGO = names(pVals_Weight_40h_Down_CC))
GOData_40h_Down_CC_DE$DEG <- GOData_40h_Down_CC_DE$Significant
GOData_40h_Down_CC_DE$Term <- Term(rownames(GOData_40h_Down_CC_DE))

ggplot(GOData_40h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_40h_Down_CC_DE,
            "topGO_weightFS_40h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 46h Up-regulated

geneList_46h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_46h_Up$Niben))
names(geneList_46h_Up) <- geneNames_Niben261

GOData_46h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_46h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_46h_Up_CC <- runTest(GOData_46h_Up_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_46h_Up_CC <- runTest(GOData_46h_Up_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_46h_Up_CC <- runTest(GOData_46h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_46h_Up_CC <- runTest(GOData_46h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_46h_Up_CC <- usedGO(GOData_46h_Up_CC)
allRes_46h_Up_CC <- GenTable(GOData_46h_Up_CC, 
                             classicFisher = resultsFis_Class_46h_Up_CC,
                             weightFisher = resultsFis_Weight_46h_Up_CC,
                             #classicKS = resultsKS_Class_46h_Up_CC,
                             #weightKS = resultsKS_Weight_46h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_46h_Up_CC))

pVals_Weight_46h_Up_CC <- score(resultsFis_Weight_46h_Up_CC)[score(resultsFis_Weight_46h_Up_CC)<=0.05]

GOData_46h_Up_CC_DE <- termStat(object = GOData_46h_Up_CC, whichGO = names(pVals_Weight_46h_Up_CC))
GOData_46h_Up_CC_DE$DEG <- GOData_46h_Up_CC_DE$Significant
GOData_46h_Up_CC_DE$Term <- Term(rownames(GOData_46h_Up_CC_DE))

ggplot(GOData_46h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_46h_Up_CC_DE,
            "topGO_weightFS_46h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 46h Down-regulated

geneList_46h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_46h_Down$Niben))
names(geneList_46h_Down) <- geneNames_Niben261

GOData_46h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_46h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_46h_Down_CC <- runTest(GOData_46h_Down_CC, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_46h_Down_CC <- runTest(GOData_46h_Down_CC, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_46h_Down_CC <- runTest(GOData_46h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_46h_Down_CC <- runTest(GOData_46h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_46h_Down_CC <- usedGO(GOData_46h_Down_CC)
allRes_46h_Down_CC <- GenTable(GOData_46h_Down_CC, 
                               classicFisher = resultsFis_Class_46h_Down_CC,
                               weightFisher = resultsFis_Weight_46h_Down_CC,
                               #classicKS = resultsKS_Class_46h_Down_CC,
                               #weightKS = resultsKS_Weight_46h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_46h_Down_CC))

pVals_Weight_46h_Down_CC <- score(resultsFis_Weight_46h_Down_CC)[score(resultsFis_Weight_46h_Down_CC)<=0.05]

GOData_46h_Down_CC_DE <- termStat(object = GOData_46h_Down_CC, whichGO = names(pVals_Weight_46h_Down_CC))
GOData_46h_Down_CC_DE$DEG <- GOData_46h_Down_CC_DE$Significant
GOData_46h_Down_CC_DE$Term <- Term(rownames(GOData_46h_Down_CC_DE))

ggplot(GOData_46h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_46h_Down_CC_DE,
            "topGO_weightFS_46h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 56h Up-regulated

geneList_56h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_56h_Up$Niben))
names(geneList_56h_Up) <- geneNames_Niben261

GOData_56h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_56h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_56h_Up_CC <- runTest(GOData_56h_Up_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_56h_Up_CC <- runTest(GOData_56h_Up_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_56h_Up_CC <- runTest(GOData_56h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_56h_Up_CC <- runTest(GOData_56h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_56h_Up_CC <- usedGO(GOData_56h_Up_CC)
allRes_56h_Up_CC <- GenTable(GOData_56h_Up_CC, 
                             classicFisher = resultsFis_Class_56h_Up_CC,
                             weightFisher = resultsFis_Weight_56h_Up_CC,
                             #classicKS = resultsKS_Class_56h_Up_CC,
                             #weightKS = resultsKS_Weight_56h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_56h_Up_CC))

pVals_Weight_56h_Up_CC <- score(resultsFis_Weight_56h_Up_CC)[score(resultsFis_Weight_56h_Up_CC)<=0.05]

GOData_56h_Up_CC_DE <- termStat(object = GOData_56h_Up_CC, whichGO = names(pVals_Weight_56h_Up_CC))
GOData_56h_Up_CC_DE$DEG <- GOData_56h_Up_CC_DE$Significant
GOData_56h_Up_CC_DE$Term <- Term(rownames(GOData_56h_Up_CC_DE))

ggplot(GOData_56h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_56h_Up_CC_DE,
            "topGO_weightFS_56h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 56h Down-regulated

geneList_56h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_56h_Down$Niben))
names(geneList_56h_Down) <- geneNames_Niben261

GOData_56h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_56h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_56h_Down_CC <- runTest(GOData_56h_Down_CC, algorithm = "weight01", 
#                                        statistic = "ks")
#resultsKS_Class_56h_Down_CC <- runTest(GOData_56h_Down_CC, algorithm = "classic", 
#                                       statistic = "ks")
resultsFis_Class_56h_Down_CC <- runTest(GOData_56h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_56h_Down_CC <- runTest(GOData_56h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_56h_Down_CC <- usedGO(GOData_56h_Down_CC)
allRes_56h_Down_CC <- GenTable(GOData_56h_Down_CC, 
                               classicFisher = resultsFis_Class_56h_Down_CC,
                               weightFisher = resultsFis_Weight_56h_Down_CC,
                               #classicKS = resultsKS_Class_56h_Down_CC,
                               #weightKS = resultsKS_Weight_56h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_56h_Down_CC))

pVals_Weight_56h_Down_CC <- score(resultsFis_Weight_56h_Down_CC)[score(resultsFis_Weight_56h_Down_CC)<=0.05]

GOData_56h_Down_CC_DE <- termStat(object = GOData_56h_Down_CC, whichGO = names(pVals_Weight_56h_Down_CC))
GOData_56h_Down_CC_DE$DEG <- GOData_56h_Down_CC_DE$Significant
GOData_56h_Down_CC_DE$Term <- Term(rownames(GOData_56h_Down_CC_DE))

ggplot(GOData_56h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_56h_Down_CC_DE,
            "topGO_weightFS_56h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

################################################################################################################
## UpSetR diagrams of CC in Cornell TimePoints


Cornell_TimePoints_Up_CC <- list("22hUp"= c(GOData_22h_Up_CC_DE$Term),
                                 "25hUp"= c(GOData_25h_Up_CC_DE$Term),
                                 "28hUp"= c(GOData_28h_Up_CC_DE$Term),
                                 "34hUp"= c(GOData_34h_Up_CC_DE$Term),
                                 "37hUp"= c(GOData_37h_Up_CC_DE$Term),
                                 "40hUp"= c(GOData_40h_Up_CC_DE$Term),
                                 "46hUp"= c(GOData_46h_Up_CC_DE$Term),
                                 "56hUp"= c(GOData_56h_Up_CC_DE$Term))

upset(fromList(Cornell_TimePoints_Up_CC), sets= rev(c("22hUp", "25hUp", "28hUp", "34hUp", "37hUp", "40hUp", "46hUp", "56hUp")), order.by = "freq", keep.order = TRUE)

Cornell_TimePoints_Down_CC <- list("22hDown"= c(GOData_22h_Down_CC_DE$Term),
                                   "25hDown"= c(GOData_25h_Down_CC_DE$Term),
                                   "28hDown"= c(GOData_28h_Down_CC_DE$Term),
                                   "34hDown"= c(GOData_34h_Down_CC_DE$Term),
                                   "37hDown"= c(GOData_37h_Down_CC_DE$Term),
                                   "40hDown"= c(GOData_40h_Down_CC_DE$Term),
                                   "46hDown"= c(GOData_46h_Down_CC_DE$Term),
                                   "56hDown"= c(GOData_56h_Down_CC_DE$Term))

upset(fromList(Cornell_TimePoints_Down_CC), sets= rev(c("22hDown", "25hDown", "28hDown", "34hDown", "37hDown", "40hDown", "46hDown", "56hDown")), order.by = "freq", keep.order = TRUE)


##################################################################################################
################################   PCA GROUPS   ##################################################

#############################  BIOLOGICAL PROCESS  ###############################################
#################################################################################################
## AA Up-regulated

geneList_AA_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_AA_Up$Niben))
names(geneList_AA_Up) <- geneNames_Niben261

GOData_AA_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_AA_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_AA_Up_BP <- runTest(GOData_AA_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_AA_Up_BP <- runTest(GOData_AA_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_AA_Up_BP <- runTest(GOData_AA_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_AA_Up_BP <- runTest(GOData_AA_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_AA_Up_BP <- usedGO(GOData_AA_Up_BP)
allRes_AA_Up_BP <- GenTable(GOData_AA_Up_BP, 
                             classicFisher = resultsFis_Class_AA_Up_BP,
                             weightFisher = resultsFis_Weight_AA_Up_BP,
                             #classicKS = resultsKS_Class_AA_Up_BP,
                             #weightKS = resultsKS_Weight_AA_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_AA_Up_BP))

pVals_Weight_AA_Up_BP <- score(resultsFis_Weight_AA_Up_BP)[score(resultsFis_Weight_AA_Up_BP)<=0.05]

GOData_AA_Up_BP_DE <- termStat(object = GOData_AA_Up_BP, whichGO = names(pVals_Weight_AA_Up_BP))
GOData_AA_Up_BP_DE$DEG <- GOData_AA_Up_BP_DE$Significant
GOData_AA_Up_BP_DE$Term <- Term(rownames(GOData_AA_Up_BP_DE))

ggplot(GOData_AA_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_AA_Up_BP_DE,
            "topGO_weightFS_AA_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## AA Down-regulated

geneList_AA_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_AA_Down$Niben))
names(geneList_AA_Down) <- geneNames_Niben261

GOData_AA_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_AA_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_AA_Down_BP <- runTest(GOData_AA_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_AA_Down_BP <- runTest(GOData_AA_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_AA_Down_BP <- runTest(GOData_AA_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_AA_Down_BP <- runTest(GOData_AA_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_AA_Down_BP <- usedGO(GOData_AA_Down_BP)
allRes_AA_Down_BP <- GenTable(GOData_AA_Down_BP, 
                               classicFisher = resultsFis_Class_AA_Down_BP,
                               weightFisher = resultsFis_Weight_AA_Down_BP,
                               #classicKS = resultsKS_Class_AA_Down_BP,
                               #weightKS = resultsKS_Weight_AA_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_AA_Down_BP))

pVals_Weight_AA_Down_BP <- score(resultsFis_Weight_AA_Down_BP)[score(resultsFis_Weight_AA_Down_BP)<=0.05]

GOData_AA_Down_BP_DE <- termStat(object = GOData_AA_Down_BP, whichGO = names(pVals_Weight_AA_Down_BP))
GOData_AA_Down_BP_DE$DEG <- GOData_AA_Down_BP_DE$Significant
GOData_AA_Down_BP_DE$Term <- Term(rownames(GOData_AA_Down_BP_DE))

ggplot(GOData_AA_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_AA_Down_BP_DE,
            "topGO_weightFS_AA_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## BB Up-regulated

geneList_BB_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_BB_Up$Niben))
names(geneList_BB_Up) <- geneNames_Niben261

GOData_BB_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_BB_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_BB_Up_BP <- runTest(GOData_BB_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_BB_Up_BP <- runTest(GOData_BB_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_BB_Up_BP <- runTest(GOData_BB_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_BB_Up_BP <- runTest(GOData_BB_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_BB_Up_BP <- usedGO(GOData_BB_Up_BP)
allRes_BB_Up_BP <- GenTable(GOData_BB_Up_BP, 
                             classicFisher = resultsFis_Class_BB_Up_BP,
                             weightFisher = resultsFis_Weight_BB_Up_BP,
                             #classicKS = resultsKS_Class_BB_Up_BP,
                             #weightKS = resultsKS_Weight_BB_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_BB_Up_BP))

pVals_Weight_BB_Up_BP <- score(resultsFis_Weight_BB_Up_BP)[score(resultsFis_Weight_BB_Up_BP)<=0.05]

GOData_BB_Up_BP_DE <- termStat(object = GOData_BB_Up_BP, whichGO = names(pVals_Weight_BB_Up_BP))
GOData_BB_Up_BP_DE$DEG <- GOData_BB_Up_BP_DE$Significant
GOData_BB_Up_BP_DE$Term <- Term(rownames(GOData_BB_Up_BP_DE))

ggplot(GOData_BB_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_BB_Up_BP_DE,
            "topGO_weightFS_BB_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## BB Down-regulated

geneList_BB_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_BB_Down$Niben))
names(geneList_BB_Down) <- geneNames_Niben261

GOData_BB_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_BB_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_BB_Down_BP <- runTest(GOData_BB_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_BB_Down_BP <- runTest(GOData_BB_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_BB_Down_BP <- runTest(GOData_BB_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_BB_Down_BP <- runTest(GOData_BB_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_BB_Down_BP <- usedGO(GOData_BB_Down_BP)
allRes_BB_Down_BP <- GenTable(GOData_BB_Down_BP, 
                               classicFisher = resultsFis_Class_BB_Down_BP,
                               weightFisher = resultsFis_Weight_BB_Down_BP,
                               #classicKS = resultsKS_Class_BB_Down_BP,
                               #weightKS = resultsKS_Weight_BB_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_BB_Down_BP))

pVals_Weight_BB_Down_BP <- score(resultsFis_Weight_BB_Down_BP)[score(resultsFis_Weight_BB_Down_BP)<=0.05]

GOData_BB_Down_BP_DE <- termStat(object = GOData_BB_Down_BP, whichGO = names(pVals_Weight_BB_Down_BP))
GOData_BB_Down_BP_DE$DEG <- GOData_BB_Down_BP_DE$Significant
GOData_BB_Down_BP_DE$Term <- Term(rownames(GOData_BB_Down_BP_DE))

ggplot(GOData_BB_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_BB_Down_BP_DE,
            "topGO_weightFS_BB_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## CC Up-regulated

geneList_CC_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_CC_Up$Niben))
names(geneList_CC_Up) <- geneNames_Niben261

GOData_CC_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_CC_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_CC_Up_BP <- runTest(GOData_CC_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_CC_Up_BP <- runTest(GOData_CC_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_CC_Up_BP <- runTest(GOData_CC_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_CC_Up_BP <- runTest(GOData_CC_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_CC_Up_BP <- usedGO(GOData_CC_Up_BP)
allRes_CC_Up_BP <- GenTable(GOData_CC_Up_BP, 
                             classicFisher = resultsFis_Class_CC_Up_BP,
                             weightFisher = resultsFis_Weight_CC_Up_BP,
                             #classicKS = resultsKS_Class_CC_Up_BP,
                             #weightKS = resultsKS_Weight_CC_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_CC_Up_BP))

pVals_Weight_CC_Up_BP <- score(resultsFis_Weight_CC_Up_BP)[score(resultsFis_Weight_CC_Up_BP)<=0.05]

GOData_CC_Up_BP_DE <- termStat(object = GOData_CC_Up_BP, whichGO = names(pVals_Weight_CC_Up_BP))
GOData_CC_Up_BP_DE$DEG <- GOData_CC_Up_BP_DE$Significant
GOData_CC_Up_BP_DE$Term <- Term(rownames(GOData_CC_Up_BP_DE))

ggplot(GOData_CC_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_CC_Up_BP_DE,
            "topGO_weightFS_CC_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## CC Down-regulated

geneList_CC_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_CC_Down$Niben))
names(geneList_CC_Down) <- geneNames_Niben261

GOData_CC_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_CC_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_CC_Down_BP <- runTest(GOData_CC_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_CC_Down_BP <- runTest(GOData_CC_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_CC_Down_BP <- runTest(GOData_CC_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_CC_Down_BP <- runTest(GOData_CC_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_CC_Down_BP <- usedGO(GOData_CC_Down_BP)
allRes_CC_Down_BP <- GenTable(GOData_CC_Down_BP, 
                               classicFisher = resultsFis_Class_CC_Down_BP,
                               weightFisher = resultsFis_Weight_CC_Down_BP,
                               #classicKS = resultsKS_Class_CC_Down_BP,
                               #weightKS = resultsKS_Weight_CC_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_CC_Down_BP))

pVals_Weight_CC_Down_BP <- score(resultsFis_Weight_CC_Down_BP)[score(resultsFis_Weight_CC_Down_BP)<=0.05]

GOData_CC_Down_BP_DE <- termStat(object = GOData_CC_Down_BP, whichGO = names(pVals_Weight_CC_Down_BP))
GOData_CC_Down_BP_DE$DEG <- GOData_CC_Down_BP_DE$Significant
GOData_CC_Down_BP_DE$Term <- Term(rownames(GOData_CC_Down_BP_DE))

ggplot(GOData_CC_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_CC_Down_BP_DE,
            "topGO_weightFS_CC_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## DD Up-regulated

geneList_DD_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_DD_Up$Niben))
names(geneList_DD_Up) <- geneNames_Niben261

GOData_DD_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_DD_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_DD_Up_BP <- runTest(GOData_DD_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_DD_Up_BP <- runTest(GOData_DD_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_DD_Up_BP <- runTest(GOData_DD_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_DD_Up_BP <- runTest(GOData_DD_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_DD_Up_BP <- usedGO(GOData_DD_Up_BP)
allRes_DD_Up_BP <- GenTable(GOData_DD_Up_BP, 
                             classicFisher = resultsFis_Class_DD_Up_BP,
                             weightFisher = resultsFis_Weight_DD_Up_BP,
                             #classicKS = resultsKS_Class_DD_Up_BP,
                             #weightKS = resultsKS_Weight_DD_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_DD_Up_BP))

pVals_Weight_DD_Up_BP <- score(resultsFis_Weight_DD_Up_BP)[score(resultsFis_Weight_DD_Up_BP)<=0.05]

GOData_DD_Up_BP_DE <- termStat(object = GOData_DD_Up_BP, whichGO = names(pVals_Weight_DD_Up_BP))
GOData_DD_Up_BP_DE$DEG <- GOData_DD_Up_BP_DE$Significant
GOData_DD_Up_BP_DE$Term <- Term(rownames(GOData_DD_Up_BP_DE))

ggplot(GOData_DD_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_DD_Up_BP_DE,
            "topGO_weightFS_DD_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## DD Down-regulated

geneList_DD_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_DD_Down$Niben))
names(geneList_DD_Down) <- geneNames_Niben261

GOData_DD_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_DD_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_DD_Down_BP <- runTest(GOData_DD_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_DD_Down_BP <- runTest(GOData_DD_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_DD_Down_BP <- runTest(GOData_DD_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_DD_Down_BP <- runTest(GOData_DD_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_DD_Down_BP <- usedGO(GOData_DD_Down_BP)
allRes_DD_Down_BP <- GenTable(GOData_DD_Down_BP, 
                               classicFisher = resultsFis_Class_DD_Down_BP,
                               weightFisher = resultsFis_Weight_DD_Down_BP,
                               #classicKS = resultsKS_Class_DD_Down_BP,
                               #weightKS = resultsKS_Weight_DD_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_DD_Down_BP))

pVals_Weight_DD_Down_BP <- score(resultsFis_Weight_DD_Down_BP)[score(resultsFis_Weight_DD_Down_BP)<=0.05]

GOData_DD_Down_BP_DE <- termStat(object = GOData_DD_Down_BP, whichGO = names(pVals_Weight_DD_Down_BP))
GOData_DD_Down_BP_DE$DEG <- GOData_DD_Down_BP_DE$Significant
GOData_DD_Down_BP_DE$Term <- Term(rownames(GOData_DD_Down_BP_DE))

ggplot(GOData_DD_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_DD_Down_BP_DE,
            "topGO_weightFS_DD_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## EE Up-regulated

geneList_EE_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_EE_Up$Niben))
names(geneList_EE_Up) <- geneNames_Niben261

GOData_EE_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_EE_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_EE_Up_BP <- runTest(GOData_EE_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_EE_Up_BP <- runTest(GOData_EE_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_EE_Up_BP <- runTest(GOData_EE_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_EE_Up_BP <- runTest(GOData_EE_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_EE_Up_BP <- usedGO(GOData_EE_Up_BP)
allRes_EE_Up_BP <- GenTable(GOData_EE_Up_BP, 
                             classicFisher = resultsFis_Class_EE_Up_BP,
                             weightFisher = resultsFis_Weight_EE_Up_BP,
                             #classicKS = resultsKS_Class_EE_Up_BP,
                             #weightKS = resultsKS_Weight_EE_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_EE_Up_BP))

pVals_Weight_EE_Up_BP <- score(resultsFis_Weight_EE_Up_BP)[score(resultsFis_Weight_EE_Up_BP)<=0.05]

GOData_EE_Up_BP_DE <- termStat(object = GOData_EE_Up_BP, whichGO = names(pVals_Weight_EE_Up_BP))
GOData_EE_Up_BP_DE$DEG <- GOData_EE_Up_BP_DE$Significant
GOData_EE_Up_BP_DE$Term <- Term(rownames(GOData_EE_Up_BP_DE))

ggplot(GOData_EE_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_EE_Up_BP_DE,
            "topGO_weightFS_EE_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## EE Down-regulated

geneList_EE_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_EE_Down$Niben))
names(geneList_EE_Down) <- geneNames_Niben261

GOData_EE_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_EE_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_EE_Down_BP <- runTest(GOData_EE_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_EE_Down_BP <- runTest(GOData_EE_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_EE_Down_BP <- runTest(GOData_EE_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_EE_Down_BP <- runTest(GOData_EE_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_EE_Down_BP <- usedGO(GOData_EE_Down_BP)
allRes_EE_Down_BP <- GenTable(GOData_EE_Down_BP, 
                               classicFisher = resultsFis_Class_EE_Down_BP,
                               weightFisher = resultsFis_Weight_EE_Down_BP,
                               #classicKS = resultsKS_Class_EE_Down_BP,
                               #weightKS = resultsKS_Weight_EE_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_EE_Down_BP))

pVals_Weight_EE_Down_BP <- score(resultsFis_Weight_EE_Down_BP)[score(resultsFis_Weight_EE_Down_BP)<=0.05]

GOData_EE_Down_BP_DE <- termStat(object = GOData_EE_Down_BP, whichGO = names(pVals_Weight_EE_Down_BP))
GOData_EE_Down_BP_DE$DEG <- GOData_EE_Down_BP_DE$Significant
GOData_EE_Down_BP_DE$Term <- Term(rownames(GOData_EE_Down_BP_DE))

ggplot(GOData_EE_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_EE_Down_BP_DE,
            "topGO_weightFS_EE_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## YY Up-regulated

geneList_YY_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_YY_Up$Niben))
names(geneList_YY_Up) <- geneNames_Niben261

GOData_YY_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_YY_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_YY_Up_BP <- runTest(GOData_YY_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_YY_Up_BP <- runTest(GOData_YY_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_YY_Up_BP <- runTest(GOData_YY_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_YY_Up_BP <- runTest(GOData_YY_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_YY_Up_BP <- usedGO(GOData_YY_Up_BP)
allRes_YY_Up_BP <- GenTable(GOData_YY_Up_BP, 
                             classicFisher = resultsFis_Class_YY_Up_BP,
                             weightFisher = resultsFis_Weight_YY_Up_BP,
                             #classicKS = resultsKS_Class_YY_Up_BP,
                             #weightKS = resultsKS_Weight_YY_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_YY_Up_BP))

pVals_Weight_YY_Up_BP <- score(resultsFis_Weight_YY_Up_BP)[score(resultsFis_Weight_YY_Up_BP)<=0.05]

GOData_YY_Up_BP_DE <- termStat(object = GOData_YY_Up_BP, whichGO = names(pVals_Weight_YY_Up_BP))
GOData_YY_Up_BP_DE$DEG <- GOData_YY_Up_BP_DE$Significant
GOData_YY_Up_BP_DE$Term <- Term(rownames(GOData_YY_Up_BP_DE))

ggplot(GOData_YY_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_YY_Up_BP_DE,
            "topGO_weightFS_YY_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## YY Down-regulated

geneList_YY_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_YY_Down$Niben))
names(geneList_YY_Down) <- geneNames_Niben261

GOData_YY_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_YY_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_YY_Down_BP <- runTest(GOData_YY_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_YY_Down_BP <- runTest(GOData_YY_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_YY_Down_BP <- runTest(GOData_YY_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_YY_Down_BP <- runTest(GOData_YY_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_YY_Down_BP <- usedGO(GOData_YY_Down_BP)
allRes_YY_Down_BP <- GenTable(GOData_YY_Down_BP, 
                               classicFisher = resultsFis_Class_YY_Down_BP,
                               weightFisher = resultsFis_Weight_YY_Down_BP,
                               #classicKS = resultsKS_Class_YY_Down_BP,
                               #weightKS = resultsKS_Weight_YY_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_YY_Down_BP))

pVals_Weight_YY_Down_BP <- score(resultsFis_Weight_YY_Down_BP)[score(resultsFis_Weight_YY_Down_BP)<=0.05]

GOData_YY_Down_BP_DE <- termStat(object = GOData_YY_Down_BP, whichGO = names(pVals_Weight_YY_Down_BP))
GOData_YY_Down_BP_DE$DEG <- GOData_YY_Down_BP_DE$Significant
GOData_YY_Down_BP_DE$Term <- Term(rownames(GOData_YY_Down_BP_DE))

ggplot(GOData_YY_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_YY_Down_BP_DE,
            "topGO_weightFS_YY_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## ZZ Up-regulated

geneList_ZZ_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_ZZ_Up$Niben))
names(geneList_ZZ_Up) <- geneNames_Niben261

GOData_ZZ_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_ZZ_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_ZZ_Up_BP <- runTest(GOData_ZZ_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_ZZ_Up_BP <- runTest(GOData_ZZ_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_ZZ_Up_BP <- runTest(GOData_ZZ_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_ZZ_Up_BP <- runTest(GOData_ZZ_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_ZZ_Up_BP <- usedGO(GOData_ZZ_Up_BP)
allRes_ZZ_Up_BP <- GenTable(GOData_ZZ_Up_BP, 
                             classicFisher = resultsFis_Class_ZZ_Up_BP,
                             weightFisher = resultsFis_Weight_ZZ_Up_BP,
                             #classicKS = resultsKS_Class_ZZ_Up_BP,
                             #weightKS = resultsKS_Weight_ZZ_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_ZZ_Up_BP))

pVals_Weight_ZZ_Up_BP <- score(resultsFis_Weight_ZZ_Up_BP)[score(resultsFis_Weight_ZZ_Up_BP)<=0.05]

GOData_ZZ_Up_BP_DE <- termStat(object = GOData_ZZ_Up_BP, whichGO = names(pVals_Weight_ZZ_Up_BP))
GOData_ZZ_Up_BP_DE$DEG <- GOData_ZZ_Up_BP_DE$Significant
GOData_ZZ_Up_BP_DE$Term <- Term(rownames(GOData_ZZ_Up_BP_DE))

ggplot(GOData_ZZ_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_ZZ_Up_BP_DE,
            "topGO_weightFS_ZZ_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## ZZ Down-regulated

geneList_ZZ_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_ZZ_Down$Niben))
names(geneList_ZZ_Down) <- geneNames_Niben261

GOData_ZZ_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_ZZ_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_ZZ_Down_BP <- runTest(GOData_ZZ_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_ZZ_Down_BP <- runTest(GOData_ZZ_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_ZZ_Down_BP <- runTest(GOData_ZZ_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_ZZ_Down_BP <- runTest(GOData_ZZ_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_ZZ_Down_BP <- usedGO(GOData_ZZ_Down_BP)
allRes_ZZ_Down_BP <- GenTable(GOData_ZZ_Down_BP, 
                               classicFisher = resultsFis_Class_ZZ_Down_BP,
                               weightFisher = resultsFis_Weight_ZZ_Down_BP,
                               #classicKS = resultsKS_Class_ZZ_Down_BP,
                               #weightKS = resultsKS_Weight_ZZ_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_ZZ_Down_BP))

pVals_Weight_ZZ_Down_BP <- score(resultsFis_Weight_ZZ_Down_BP)[score(resultsFis_Weight_ZZ_Down_BP)<=0.05]

GOData_ZZ_Down_BP_DE <- termStat(object = GOData_ZZ_Down_BP, whichGO = names(pVals_Weight_ZZ_Down_BP))
GOData_ZZ_Down_BP_DE$DEG <- GOData_ZZ_Down_BP_DE$Significant
GOData_ZZ_Down_BP_DE$Term <- Term(rownames(GOData_ZZ_Down_BP_DE))

ggplot(GOData_ZZ_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_ZZ_Down_BP_DE,
            "topGO_weightFS_ZZ_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

################################################################################################################
## UpSetR diagrams of BP in Cornell TimePoints


Cornell_PCA_Up_BP <- list("AAUp"= c(GOData_AA_Up_BP_DE$Term),
                          "BBUp"= c(GOData_BB_Up_BP_DE$Term),
                          "YYUp"= c(GOData_YY_Up_BP_DE$Term),
                          "ZZUp"= c(GOData_ZZ_Up_BP_DE$Term),
                          "CCUp"= c(GOData_CC_Up_BP_DE$Term),
                          "DDUp"= c(GOData_DD_Up_BP_DE$Term),
                          "EEUp"= c(GOData_EE_Up_BP_DE$Term),
                          "56hUp"= c(GOData_56h_Up_BP_DE$Term))

upset(fromList(Cornell_PCA_Up_BP), sets= rev(c("AAUp", "BBUp", "YYUp", "ZZUp", "CCUp", "DDUp", "EEUp", "56hUp")), order.by = "freq", keep.order = TRUE)

Cornell_PCA_Down_BP <- list("AADown"= c(GOData_AA_Down_BP_DE$Term),
                          "BBDown"= c(GOData_BB_Down_BP_DE$Term),
                          "YYDown"= c(GOData_YY_Down_BP_DE$Term),
                          "ZZDown"= c(GOData_ZZ_Down_BP_DE$Term),
                          "CCDown"= c(GOData_CC_Down_BP_DE$Term),
                          "DDDown"= c(GOData_DD_Down_BP_DE$Term),
                          "EEDown"= c(GOData_EE_Down_BP_DE$Term),
                          "56hDown"= c(GOData_56h_Down_BP_DE$Term))

upset(fromList(Cornell_PCA_Down_BP), sets= rev(c("AADown", "BBDown", "YYDown", "ZZDown", "CCDown", "DDDown", "EEDown", "56hDown")), order.by = "freq", keep.order = TRUE)


#############################  MOLECULAR FUNCTION  ###############################################
#################################################################################################
## AA Up-regulated

geneList_AA_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_AA_Up$Niben))
names(geneList_AA_Up) <- geneNames_Niben261

GOData_AA_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_AA_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_AA_Up_MF <- runTest(GOData_AA_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_AA_Up_MF <- runTest(GOData_AA_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_AA_Up_MF <- runTest(GOData_AA_Up_MF, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_AA_Up_MF <- runTest(GOData_AA_Up_MF, algorithm = "weight01",
                                      statistic = "fisher")

allGO_AA_Up_MF <- usedGO(GOData_AA_Up_MF)
allRes_AA_Up_MF <- GenTable(GOData_AA_Up_MF, 
                            classicFisher = resultsFis_Class_AA_Up_MF,
                            weightFisher = resultsFis_Weight_AA_Up_MF,
                            #classicKS = resultsKS_Class_AA_Up_MF,
                            #weightKS = resultsKS_Weight_AA_Up_MF,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_AA_Up_MF))

pVals_Weight_AA_Up_MF <- score(resultsFis_Weight_AA_Up_MF)[score(resultsFis_Weight_AA_Up_MF)<=0.05]

GOData_AA_Up_MF_DE <- termStat(object = GOData_AA_Up_MF, whichGO = names(pVals_Weight_AA_Up_MF))
GOData_AA_Up_MF_DE$DEG <- GOData_AA_Up_MF_DE$Significant
GOData_AA_Up_MF_DE$Term <- Term(rownames(GOData_AA_Up_MF_DE))

ggplot(GOData_AA_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_AA_Up_MF_DE,
            "topGO_weightFS_AA_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## AA Down-regulated

geneList_AA_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_AA_Down$Niben))
names(geneList_AA_Down) <- geneNames_Niben261

GOData_AA_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_AA_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_AA_Down_MF <- runTest(GOData_AA_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_AA_Down_MF <- runTest(GOData_AA_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_AA_Down_MF <- runTest(GOData_AA_Down_MF, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_AA_Down_MF <- runTest(GOData_AA_Down_MF, algorithm = "weight01",
                                        statistic = "fisher")

allGO_AA_Down_MF <- usedGO(GOData_AA_Down_MF)
allRes_AA_Down_MF <- GenTable(GOData_AA_Down_MF, 
                              classicFisher = resultsFis_Class_AA_Down_MF,
                              weightFisher = resultsFis_Weight_AA_Down_MF,
                              #classicKS = resultsKS_Class_AA_Down_MF,
                              #weightKS = resultsKS_Weight_AA_Down_MF,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_AA_Down_MF))

pVals_Weight_AA_Down_MF <- score(resultsFis_Weight_AA_Down_MF)[score(resultsFis_Weight_AA_Down_MF)<=0.05]

GOData_AA_Down_MF_DE <- termStat(object = GOData_AA_Down_MF, whichGO = names(pVals_Weight_AA_Down_MF))
GOData_AA_Down_MF_DE$DEG <- GOData_AA_Down_MF_DE$Significant
GOData_AA_Down_MF_DE$Term <- Term(rownames(GOData_AA_Down_MF_DE))

ggplot(GOData_AA_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_AA_Down_MF_DE,
            "topGO_weightFS_AA_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## BB Up-regulated

geneList_BB_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_BB_Up$Niben))
names(geneList_BB_Up) <- geneNames_Niben261

GOData_BB_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_BB_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_BB_Up_MF <- runTest(GOData_BB_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_BB_Up_MF <- runTest(GOData_BB_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_BB_Up_MF <- runTest(GOData_BB_Up_MF, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_BB_Up_MF <- runTest(GOData_BB_Up_MF, algorithm = "weight01",
                                      statistic = "fisher")

allGO_BB_Up_MF <- usedGO(GOData_BB_Up_MF)
allRes_BB_Up_MF <- GenTable(GOData_BB_Up_MF, 
                            classicFisher = resultsFis_Class_BB_Up_MF,
                            weightFisher = resultsFis_Weight_BB_Up_MF,
                            #classicKS = resultsKS_Class_BB_Up_MF,
                            #weightKS = resultsKS_Weight_BB_Up_MF,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_BB_Up_MF))

pVals_Weight_BB_Up_MF <- score(resultsFis_Weight_BB_Up_MF)[score(resultsFis_Weight_BB_Up_MF)<=0.05]

GOData_BB_Up_MF_DE <- termStat(object = GOData_BB_Up_MF, whichGO = names(pVals_Weight_BB_Up_MF))
GOData_BB_Up_MF_DE$DEG <- GOData_BB_Up_MF_DE$Significant
GOData_BB_Up_MF_DE$Term <- Term(rownames(GOData_BB_Up_MF_DE))

ggplot(GOData_BB_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_BB_Up_MF_DE,
            "topGO_weightFS_BB_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## BB Down-regulated

geneList_BB_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_BB_Down$Niben))
names(geneList_BB_Down) <- geneNames_Niben261

GOData_BB_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_BB_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_BB_Down_MF <- runTest(GOData_BB_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_BB_Down_MF <- runTest(GOData_BB_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_BB_Down_MF <- runTest(GOData_BB_Down_MF, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_BB_Down_MF <- runTest(GOData_BB_Down_MF, algorithm = "weight01",
                                        statistic = "fisher")

allGO_BB_Down_MF <- usedGO(GOData_BB_Down_MF)
allRes_BB_Down_MF <- GenTable(GOData_BB_Down_MF, 
                              classicFisher = resultsFis_Class_BB_Down_MF,
                              weightFisher = resultsFis_Weight_BB_Down_MF,
                              #classicKS = resultsKS_Class_BB_Down_MF,
                              #weightKS = resultsKS_Weight_BB_Down_MF,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_BB_Down_MF))

pVals_Weight_BB_Down_MF <- score(resultsFis_Weight_BB_Down_MF)[score(resultsFis_Weight_BB_Down_MF)<=0.05]

GOData_BB_Down_MF_DE <- termStat(object = GOData_BB_Down_MF, whichGO = names(pVals_Weight_BB_Down_MF))
GOData_BB_Down_MF_DE$DEG <- GOData_BB_Down_MF_DE$Significant
GOData_BB_Down_MF_DE$Term <- Term(rownames(GOData_BB_Down_MF_DE))

ggplot(GOData_BB_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_BB_Down_MF_DE,
            "topGO_weightFS_BB_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## CC Up-regulated

geneList_CC_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_CC_Up$Niben))
names(geneList_CC_Up) <- geneNames_Niben261

GOData_CC_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_CC_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_CC_Up_MF <- runTest(GOData_CC_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_CC_Up_MF <- runTest(GOData_CC_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_CC_Up_MF <- runTest(GOData_CC_Up_MF, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_CC_Up_MF <- runTest(GOData_CC_Up_MF, algorithm = "weight01",
                                      statistic = "fisher")

allGO_CC_Up_MF <- usedGO(GOData_CC_Up_MF)
allRes_CC_Up_MF <- GenTable(GOData_CC_Up_MF, 
                            classicFisher = resultsFis_Class_CC_Up_MF,
                            weightFisher = resultsFis_Weight_CC_Up_MF,
                            #classicKS = resultsKS_Class_CC_Up_MF,
                            #weightKS = resultsKS_Weight_CC_Up_MF,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_CC_Up_MF))

pVals_Weight_CC_Up_MF <- score(resultsFis_Weight_CC_Up_MF)[score(resultsFis_Weight_CC_Up_MF)<=0.05]

GOData_CC_Up_MF_DE <- termStat(object = GOData_CC_Up_MF, whichGO = names(pVals_Weight_CC_Up_MF))
GOData_CC_Up_MF_DE$DEG <- GOData_CC_Up_MF_DE$Significant
GOData_CC_Up_MF_DE$Term <- Term(rownames(GOData_CC_Up_MF_DE))

ggplot(GOData_CC_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_CC_Up_MF_DE,
            "topGO_weightFS_CC_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## CC Down-regulated

geneList_CC_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_CC_Down$Niben))
names(geneList_CC_Down) <- geneNames_Niben261

GOData_CC_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_CC_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_CC_Down_MF <- runTest(GOData_CC_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_CC_Down_MF <- runTest(GOData_CC_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_CC_Down_MF <- runTest(GOData_CC_Down_MF, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_CC_Down_MF <- runTest(GOData_CC_Down_MF, algorithm = "weight01",
                                        statistic = "fisher")

allGO_CC_Down_MF <- usedGO(GOData_CC_Down_MF)
allRes_CC_Down_MF <- GenTable(GOData_CC_Down_MF, 
                              classicFisher = resultsFis_Class_CC_Down_MF,
                              weightFisher = resultsFis_Weight_CC_Down_MF,
                              #classicKS = resultsKS_Class_CC_Down_MF,
                              #weightKS = resultsKS_Weight_CC_Down_MF,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_CC_Down_MF))

pVals_Weight_CC_Down_MF <- score(resultsFis_Weight_CC_Down_MF)[score(resultsFis_Weight_CC_Down_MF)<=0.05]

GOData_CC_Down_MF_DE <- termStat(object = GOData_CC_Down_MF, whichGO = names(pVals_Weight_CC_Down_MF))
GOData_CC_Down_MF_DE$DEG <- GOData_CC_Down_MF_DE$Significant
GOData_CC_Down_MF_DE$Term <- Term(rownames(GOData_CC_Down_MF_DE))

ggplot(GOData_CC_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_CC_Down_MF_DE,
            "topGO_weightFS_CC_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## DD Up-regulated

geneList_DD_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_DD_Up$Niben))
names(geneList_DD_Up) <- geneNames_Niben261

GOData_DD_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_DD_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_DD_Up_MF <- runTest(GOData_DD_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_DD_Up_MF <- runTest(GOData_DD_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_DD_Up_MF <- runTest(GOData_DD_Up_MF, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_DD_Up_MF <- runTest(GOData_DD_Up_MF, algorithm = "weight01",
                                      statistic = "fisher")

allGO_DD_Up_MF <- usedGO(GOData_DD_Up_MF)
allRes_DD_Up_MF <- GenTable(GOData_DD_Up_MF, 
                            classicFisher = resultsFis_Class_DD_Up_MF,
                            weightFisher = resultsFis_Weight_DD_Up_MF,
                            #classicKS = resultsKS_Class_DD_Up_MF,
                            #weightKS = resultsKS_Weight_DD_Up_MF,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_DD_Up_MF))

pVals_Weight_DD_Up_MF <- score(resultsFis_Weight_DD_Up_MF)[score(resultsFis_Weight_DD_Up_MF)<=0.05]

GOData_DD_Up_MF_DE <- termStat(object = GOData_DD_Up_MF, whichGO = names(pVals_Weight_DD_Up_MF))
GOData_DD_Up_MF_DE$DEG <- GOData_DD_Up_MF_DE$Significant
GOData_DD_Up_MF_DE$Term <- Term(rownames(GOData_DD_Up_MF_DE))

ggplot(GOData_DD_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_DD_Up_MF_DE,
            "topGO_weightFS_DD_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## DD Down-regulated

geneList_DD_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_DD_Down$Niben))
names(geneList_DD_Down) <- geneNames_Niben261

GOData_DD_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_DD_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_DD_Down_MF <- runTest(GOData_DD_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_DD_Down_MF <- runTest(GOData_DD_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_DD_Down_MF <- runTest(GOData_DD_Down_MF, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_DD_Down_MF <- runTest(GOData_DD_Down_MF, algorithm = "weight01",
                                        statistic = "fisher")

allGO_DD_Down_MF <- usedGO(GOData_DD_Down_MF)
allRes_DD_Down_MF <- GenTable(GOData_DD_Down_MF, 
                              classicFisher = resultsFis_Class_DD_Down_MF,
                              weightFisher = resultsFis_Weight_DD_Down_MF,
                              #classicKS = resultsKS_Class_DD_Down_MF,
                              #weightKS = resultsKS_Weight_DD_Down_MF,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_DD_Down_MF))

pVals_Weight_DD_Down_MF <- score(resultsFis_Weight_DD_Down_MF)[score(resultsFis_Weight_DD_Down_MF)<=0.05]

GOData_DD_Down_MF_DE <- termStat(object = GOData_DD_Down_MF, whichGO = names(pVals_Weight_DD_Down_MF))
GOData_DD_Down_MF_DE$DEG <- GOData_DD_Down_MF_DE$Significant
GOData_DD_Down_MF_DE$Term <- Term(rownames(GOData_DD_Down_MF_DE))

ggplot(GOData_DD_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_DD_Down_MF_DE,
            "topGO_weightFS_DD_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## EE Up-regulated

geneList_EE_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_EE_Up$Niben))
names(geneList_EE_Up) <- geneNames_Niben261

GOData_EE_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_EE_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_EE_Up_MF <- runTest(GOData_EE_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_EE_Up_MF <- runTest(GOData_EE_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_EE_Up_MF <- runTest(GOData_EE_Up_MF, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_EE_Up_MF <- runTest(GOData_EE_Up_MF, algorithm = "weight01",
                                      statistic = "fisher")

allGO_EE_Up_MF <- usedGO(GOData_EE_Up_MF)
allRes_EE_Up_MF <- GenTable(GOData_EE_Up_MF, 
                            classicFisher = resultsFis_Class_EE_Up_MF,
                            weightFisher = resultsFis_Weight_EE_Up_MF,
                            #classicKS = resultsKS_Class_EE_Up_MF,
                            #weightKS = resultsKS_Weight_EE_Up_MF,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_EE_Up_MF))

pVals_Weight_EE_Up_MF <- score(resultsFis_Weight_EE_Up_MF)[score(resultsFis_Weight_EE_Up_MF)<=0.05]

GOData_EE_Up_MF_DE <- termStat(object = GOData_EE_Up_MF, whichGO = names(pVals_Weight_EE_Up_MF))
GOData_EE_Up_MF_DE$DEG <- GOData_EE_Up_MF_DE$Significant
GOData_EE_Up_MF_DE$Term <- Term(rownames(GOData_EE_Up_MF_DE))

ggplot(GOData_EE_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_EE_Up_MF_DE,
            "topGO_weightFS_EE_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## EE Down-regulated

geneList_EE_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_EE_Down$Niben))
names(geneList_EE_Down) <- geneNames_Niben261

GOData_EE_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_EE_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_EE_Down_MF <- runTest(GOData_EE_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_EE_Down_MF <- runTest(GOData_EE_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_EE_Down_MF <- runTest(GOData_EE_Down_MF, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_EE_Down_MF <- runTest(GOData_EE_Down_MF, algorithm = "weight01",
                                        statistic = "fisher")

allGO_EE_Down_MF <- usedGO(GOData_EE_Down_MF)
allRes_EE_Down_MF <- GenTable(GOData_EE_Down_MF, 
                              classicFisher = resultsFis_Class_EE_Down_MF,
                              weightFisher = resultsFis_Weight_EE_Down_MF,
                              #classicKS = resultsKS_Class_EE_Down_MF,
                              #weightKS = resultsKS_Weight_EE_Down_MF,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_EE_Down_MF))

pVals_Weight_EE_Down_MF <- score(resultsFis_Weight_EE_Down_MF)[score(resultsFis_Weight_EE_Down_MF)<=0.05]

GOData_EE_Down_MF_DE <- termStat(object = GOData_EE_Down_MF, whichGO = names(pVals_Weight_EE_Down_MF))
GOData_EE_Down_MF_DE$DEG <- GOData_EE_Down_MF_DE$Significant
GOData_EE_Down_MF_DE$Term <- Term(rownames(GOData_EE_Down_MF_DE))

ggplot(GOData_EE_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_EE_Down_MF_DE,
            "topGO_weightFS_EE_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## YY Up-regulated

geneList_YY_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_YY_Up$Niben))
names(geneList_YY_Up) <- geneNames_Niben261

GOData_YY_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_YY_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_YY_Up_MF <- runTest(GOData_YY_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_YY_Up_MF <- runTest(GOData_YY_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_YY_Up_MF <- runTest(GOData_YY_Up_MF, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_YY_Up_MF <- runTest(GOData_YY_Up_MF, algorithm = "weight01",
                                      statistic = "fisher")

allGO_YY_Up_MF <- usedGO(GOData_YY_Up_MF)
allRes_YY_Up_MF <- GenTable(GOData_YY_Up_MF, 
                            classicFisher = resultsFis_Class_YY_Up_MF,
                            weightFisher = resultsFis_Weight_YY_Up_MF,
                            #classicKS = resultsKS_Class_YY_Up_MF,
                            #weightKS = resultsKS_Weight_YY_Up_MF,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_YY_Up_MF))

pVals_Weight_YY_Up_MF <- score(resultsFis_Weight_YY_Up_MF)[score(resultsFis_Weight_YY_Up_MF)<=0.05]

GOData_YY_Up_MF_DE <- termStat(object = GOData_YY_Up_MF, whichGO = names(pVals_Weight_YY_Up_MF))
GOData_YY_Up_MF_DE$DEG <- GOData_YY_Up_MF_DE$Significant
GOData_YY_Up_MF_DE$Term <- Term(rownames(GOData_YY_Up_MF_DE))

ggplot(GOData_YY_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_YY_Up_MF_DE,
            "topGO_weightFS_YY_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## YY Down-regulated

geneList_YY_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_YY_Down$Niben))
names(geneList_YY_Down) <- geneNames_Niben261

GOData_YY_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_YY_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_YY_Down_MF <- runTest(GOData_YY_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_YY_Down_MF <- runTest(GOData_YY_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_YY_Down_MF <- runTest(GOData_YY_Down_MF, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_YY_Down_MF <- runTest(GOData_YY_Down_MF, algorithm = "weight01",
                                        statistic = "fisher")

allGO_YY_Down_MF <- usedGO(GOData_YY_Down_MF)
allRes_YY_Down_MF <- GenTable(GOData_YY_Down_MF, 
                              classicFisher = resultsFis_Class_YY_Down_MF,
                              weightFisher = resultsFis_Weight_YY_Down_MF,
                              #classicKS = resultsKS_Class_YY_Down_MF,
                              #weightKS = resultsKS_Weight_YY_Down_MF,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_YY_Down_MF))

pVals_Weight_YY_Down_MF <- score(resultsFis_Weight_YY_Down_MF)[score(resultsFis_Weight_YY_Down_MF)<=0.05]

GOData_YY_Down_MF_DE <- termStat(object = GOData_YY_Down_MF, whichGO = names(pVals_Weight_YY_Down_MF))
GOData_YY_Down_MF_DE$DEG <- GOData_YY_Down_MF_DE$Significant
GOData_YY_Down_MF_DE$Term <- Term(rownames(GOData_YY_Down_MF_DE))

ggplot(GOData_YY_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_YY_Down_MF_DE,
            "topGO_weightFS_YY_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## ZZ Up-regulated

geneList_ZZ_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_ZZ_Up$Niben))
names(geneList_ZZ_Up) <- geneNames_Niben261

GOData_ZZ_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_ZZ_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_ZZ_Up_MF <- runTest(GOData_ZZ_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_ZZ_Up_MF <- runTest(GOData_ZZ_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_ZZ_Up_MF <- runTest(GOData_ZZ_Up_MF, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_ZZ_Up_MF <- runTest(GOData_ZZ_Up_MF, algorithm = "weight01",
                                      statistic = "fisher")

allGO_ZZ_Up_MF <- usedGO(GOData_ZZ_Up_MF)
allRes_ZZ_Up_MF <- GenTable(GOData_ZZ_Up_MF, 
                            classicFisher = resultsFis_Class_ZZ_Up_MF,
                            weightFisher = resultsFis_Weight_ZZ_Up_MF,
                            #classicKS = resultsKS_Class_ZZ_Up_MF,
                            #weightKS = resultsKS_Weight_ZZ_Up_MF,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_ZZ_Up_MF))

pVals_Weight_ZZ_Up_MF <- score(resultsFis_Weight_ZZ_Up_MF)[score(resultsFis_Weight_ZZ_Up_MF)<=0.05]

GOData_ZZ_Up_MF_DE <- termStat(object = GOData_ZZ_Up_MF, whichGO = names(pVals_Weight_ZZ_Up_MF))
GOData_ZZ_Up_MF_DE$DEG <- GOData_ZZ_Up_MF_DE$Significant
GOData_ZZ_Up_MF_DE$Term <- Term(rownames(GOData_ZZ_Up_MF_DE))

ggplot(GOData_ZZ_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_ZZ_Up_MF_DE,
            "topGO_weightFS_ZZ_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## ZZ Down-regulated

geneList_ZZ_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_ZZ_Down$Niben))
names(geneList_ZZ_Down) <- geneNames_Niben261

GOData_ZZ_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_ZZ_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_ZZ_Down_MF <- runTest(GOData_ZZ_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_ZZ_Down_MF <- runTest(GOData_ZZ_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_ZZ_Down_MF <- runTest(GOData_ZZ_Down_MF, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_ZZ_Down_MF <- runTest(GOData_ZZ_Down_MF, algorithm = "weight01",
                                        statistic = "fisher")

allGO_ZZ_Down_MF <- usedGO(GOData_ZZ_Down_MF)
allRes_ZZ_Down_MF <- GenTable(GOData_ZZ_Down_MF, 
                              classicFisher = resultsFis_Class_ZZ_Down_MF,
                              weightFisher = resultsFis_Weight_ZZ_Down_MF,
                              #classicKS = resultsKS_Class_ZZ_Down_MF,
                              #weightKS = resultsKS_Weight_ZZ_Down_MF,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_ZZ_Down_MF))

pVals_Weight_ZZ_Down_MF <- score(resultsFis_Weight_ZZ_Down_MF)[score(resultsFis_Weight_ZZ_Down_MF)<=0.05]

GOData_ZZ_Down_MF_DE <- termStat(object = GOData_ZZ_Down_MF, whichGO = names(pVals_Weight_ZZ_Down_MF))
GOData_ZZ_Down_MF_DE$DEG <- GOData_ZZ_Down_MF_DE$Significant
GOData_ZZ_Down_MF_DE$Term <- Term(rownames(GOData_ZZ_Down_MF_DE))

ggplot(GOData_ZZ_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_ZZ_Down_MF_DE,
            "topGO_weightFS_ZZ_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

################################################################################################################
## UpSetR diagrams of MF in Cornell TimePoints


Cornell_PCA_Up_MF <- list("AAUp"= c(GOData_AA_Up_MF_DE$Term),
                          "BBUp"= c(GOData_BB_Up_MF_DE$Term),
                          "YYUp"= c(GOData_YY_Up_MF_DE$Term),
                          "ZZUp"= c(GOData_ZZ_Up_MF_DE$Term),
                          "CCUp"= c(GOData_CC_Up_MF_DE$Term),
                          "DDUp"= c(GOData_DD_Up_MF_DE$Term),
                          "EEUp"= c(GOData_EE_Up_MF_DE$Term),
                          "56hUp"= c(GOData_56h_Up_MF_DE$Term))

upset(fromList(Cornell_PCA_Up_MF), sets= rev(c("AAUp", "BBUp", "YYUp", "ZZUp", "CCUp", "DDUp", "EEUp", "56hUp")), order.by = "freq", keep.order = TRUE)

Cornell_PCA_Down_MF <- list("AADown"= c(GOData_AA_Down_MF_DE$Term),
                            "BBDown"= c(GOData_BB_Down_MF_DE$Term),
                            "YYDown"= c(GOData_YY_Down_MF_DE$Term),
                            "ZZDown"= c(GOData_ZZ_Down_MF_DE$Term),
                            "CCDown"= c(GOData_CC_Down_MF_DE$Term),
                            "DDDown"= c(GOData_DD_Down_MF_DE$Term),
                            "EEDown"= c(GOData_EE_Down_MF_DE$Term),
                            "56hDown"= c(GOData_56h_Down_MF_DE$Term))

upset(fromList(Cornell_PCA_Down_MF), sets= rev(c("AADown", "BBDown", "YYDown", "ZZDown", "CCDown", "DDDown", "EEDown", "56hDown")), order.by = "freq", keep.order = TRUE)



#############################  CELLULAR COMPONENT  ###############################################
#################################################################################################
## AA Up-regulated

geneList_AA_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_AA_Up$Niben))
names(geneList_AA_Up) <- geneNames_Niben261

GOData_AA_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_AA_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_AA_Up_CC <- runTest(GOData_AA_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_AA_Up_CC <- runTest(GOData_AA_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_AA_Up_CC <- runTest(GOData_AA_Up_CC, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_AA_Up_CC <- runTest(GOData_AA_Up_CC, algorithm = "weight01",
                                      statistic = "fisher")

allGO_AA_Up_CC <- usedGO(GOData_AA_Up_CC)
allRes_AA_Up_CC <- GenTable(GOData_AA_Up_CC, 
                            classicFisher = resultsFis_Class_AA_Up_CC,
                            weightFisher = resultsFis_Weight_AA_Up_CC,
                            #classicKS = resultsKS_Class_AA_Up_CC,
                            #weightKS = resultsKS_Weight_AA_Up_CC,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_AA_Up_CC))

pVals_Weight_AA_Up_CC <- score(resultsFis_Weight_AA_Up_CC)[score(resultsFis_Weight_AA_Up_CC)<=0.05]

GOData_AA_Up_CC_DE <- termStat(object = GOData_AA_Up_CC, whichGO = names(pVals_Weight_AA_Up_CC))
GOData_AA_Up_CC_DE$DEG <- GOData_AA_Up_CC_DE$Significant
GOData_AA_Up_CC_DE$Term <- Term(rownames(GOData_AA_Up_CC_DE))

ggplot(GOData_AA_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_AA_Up_CC_DE,
            "topGO_weightFS_AA_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## AA Down-regulated

geneList_AA_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_AA_Down$Niben))
names(geneList_AA_Down) <- geneNames_Niben261

GOData_AA_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_AA_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_AA_Down_CC <- runTest(GOData_AA_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_AA_Down_CC <- runTest(GOData_AA_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_AA_Down_CC <- runTest(GOData_AA_Down_CC, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_AA_Down_CC <- runTest(GOData_AA_Down_CC, algorithm = "weight01",
                                        statistic = "fisher")

allGO_AA_Down_CC <- usedGO(GOData_AA_Down_CC)
allRes_AA_Down_CC <- GenTable(GOData_AA_Down_CC, 
                              classicFisher = resultsFis_Class_AA_Down_CC,
                              weightFisher = resultsFis_Weight_AA_Down_CC,
                              #classicKS = resultsKS_Class_AA_Down_CC,
                              #weightKS = resultsKS_Weight_AA_Down_CC,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_AA_Down_CC))

pVals_Weight_AA_Down_CC <- score(resultsFis_Weight_AA_Down_CC)[score(resultsFis_Weight_AA_Down_CC)<=0.05]

GOData_AA_Down_CC_DE <- termStat(object = GOData_AA_Down_CC, whichGO = names(pVals_Weight_AA_Down_CC))
GOData_AA_Down_CC_DE$DEG <- GOData_AA_Down_CC_DE$Significant
GOData_AA_Down_CC_DE$Term <- Term(rownames(GOData_AA_Down_CC_DE))

ggplot(GOData_AA_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_AA_Down_CC_DE,
            "topGO_weightFS_AA_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## BB Up-regulated

geneList_BB_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_BB_Up$Niben))
names(geneList_BB_Up) <- geneNames_Niben261

GOData_BB_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_BB_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_BB_Up_CC <- runTest(GOData_BB_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_BB_Up_CC <- runTest(GOData_BB_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_BB_Up_CC <- runTest(GOData_BB_Up_CC, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_BB_Up_CC <- runTest(GOData_BB_Up_CC, algorithm = "weight01",
                                      statistic = "fisher")

allGO_BB_Up_CC <- usedGO(GOData_BB_Up_CC)
allRes_BB_Up_CC <- GenTable(GOData_BB_Up_CC, 
                            classicFisher = resultsFis_Class_BB_Up_CC,
                            weightFisher = resultsFis_Weight_BB_Up_CC,
                            #classicKS = resultsKS_Class_BB_Up_CC,
                            #weightKS = resultsKS_Weight_BB_Up_CC,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_BB_Up_CC))

pVals_Weight_BB_Up_CC <- score(resultsFis_Weight_BB_Up_CC)[score(resultsFis_Weight_BB_Up_CC)<=0.05]

GOData_BB_Up_CC_DE <- termStat(object = GOData_BB_Up_CC, whichGO = names(pVals_Weight_BB_Up_CC))
GOData_BB_Up_CC_DE$DEG <- GOData_BB_Up_CC_DE$Significant
GOData_BB_Up_CC_DE$Term <- Term(rownames(GOData_BB_Up_CC_DE))

ggplot(GOData_BB_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_BB_Up_CC_DE,
            "topGO_weightFS_BB_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## BB Down-regulated

geneList_BB_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_BB_Down$Niben))
names(geneList_BB_Down) <- geneNames_Niben261

GOData_BB_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_BB_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_BB_Down_CC <- runTest(GOData_BB_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_BB_Down_CC <- runTest(GOData_BB_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_BB_Down_CC <- runTest(GOData_BB_Down_CC, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_BB_Down_CC <- runTest(GOData_BB_Down_CC, algorithm = "weight01",
                                        statistic = "fisher")

allGO_BB_Down_CC <- usedGO(GOData_BB_Down_CC)
allRes_BB_Down_CC <- GenTable(GOData_BB_Down_CC, 
                              classicFisher = resultsFis_Class_BB_Down_CC,
                              weightFisher = resultsFis_Weight_BB_Down_CC,
                              #classicKS = resultsKS_Class_BB_Down_CC,
                              #weightKS = resultsKS_Weight_BB_Down_CC,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_BB_Down_CC))

pVals_Weight_BB_Down_CC <- score(resultsFis_Weight_BB_Down_CC)[score(resultsFis_Weight_BB_Down_CC)<=0.05]

GOData_BB_Down_CC_DE <- termStat(object = GOData_BB_Down_CC, whichGO = names(pVals_Weight_BB_Down_CC))
GOData_BB_Down_CC_DE$DEG <- GOData_BB_Down_CC_DE$Significant
GOData_BB_Down_CC_DE$Term <- Term(rownames(GOData_BB_Down_CC_DE))

ggplot(GOData_BB_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_BB_Down_CC_DE,
            "topGO_weightFS_BB_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## CC Up-regulated

geneList_CC_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_CC_Up$Niben))
names(geneList_CC_Up) <- geneNames_Niben261

GOData_CC_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_CC_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_CC_Up_CC <- runTest(GOData_CC_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_CC_Up_CC <- runTest(GOData_CC_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_CC_Up_CC <- runTest(GOData_CC_Up_CC, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_CC_Up_CC <- runTest(GOData_CC_Up_CC, algorithm = "weight01",
                                      statistic = "fisher")

allGO_CC_Up_CC <- usedGO(GOData_CC_Up_CC)
allRes_CC_Up_CC <- GenTable(GOData_CC_Up_CC, 
                            classicFisher = resultsFis_Class_CC_Up_CC,
                            weightFisher = resultsFis_Weight_CC_Up_CC,
                            #classicKS = resultsKS_Class_CC_Up_CC,
                            #weightKS = resultsKS_Weight_CC_Up_CC,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_CC_Up_CC))

pVals_Weight_CC_Up_CC <- score(resultsFis_Weight_CC_Up_CC)[score(resultsFis_Weight_CC_Up_CC)<=0.05]

GOData_CC_Up_CC_DE <- termStat(object = GOData_CC_Up_CC, whichGO = names(pVals_Weight_CC_Up_CC))
GOData_CC_Up_CC_DE$DEG <- GOData_CC_Up_CC_DE$Significant
GOData_CC_Up_CC_DE$Term <- Term(rownames(GOData_CC_Up_CC_DE))

ggplot(GOData_CC_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_CC_Up_CC_DE,
            "topGO_weightFS_CC_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## CC Down-regulated

geneList_CC_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_CC_Down$Niben))
names(geneList_CC_Down) <- geneNames_Niben261

GOData_CC_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_CC_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_CC_Down_CC <- runTest(GOData_CC_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_CC_Down_CC <- runTest(GOData_CC_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_CC_Down_CC <- runTest(GOData_CC_Down_CC, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_CC_Down_CC <- runTest(GOData_CC_Down_CC, algorithm = "weight01",
                                        statistic = "fisher")

allGO_CC_Down_CC <- usedGO(GOData_CC_Down_CC)
allRes_CC_Down_CC <- GenTable(GOData_CC_Down_CC, 
                              classicFisher = resultsFis_Class_CC_Down_CC,
                              weightFisher = resultsFis_Weight_CC_Down_CC,
                              #classicKS = resultsKS_Class_CC_Down_CC,
                              #weightKS = resultsKS_Weight_CC_Down_CC,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_CC_Down_CC))

pVals_Weight_CC_Down_CC <- score(resultsFis_Weight_CC_Down_CC)[score(resultsFis_Weight_CC_Down_CC)<=0.05]

GOData_CC_Down_CC_DE <- termStat(object = GOData_CC_Down_CC, whichGO = names(pVals_Weight_CC_Down_CC))
GOData_CC_Down_CC_DE$DEG <- GOData_CC_Down_CC_DE$Significant
GOData_CC_Down_CC_DE$Term <- Term(rownames(GOData_CC_Down_CC_DE))

ggplot(GOData_CC_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_CC_Down_CC_DE,
            "topGO_weightFS_CC_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## DD Up-regulated

geneList_DD_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_DD_Up$Niben))
names(geneList_DD_Up) <- geneNames_Niben261

GOData_DD_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_DD_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_DD_Up_CC <- runTest(GOData_DD_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_DD_Up_CC <- runTest(GOData_DD_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_DD_Up_CC <- runTest(GOData_DD_Up_CC, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_DD_Up_CC <- runTest(GOData_DD_Up_CC, algorithm = "weight01",
                                      statistic = "fisher")

allGO_DD_Up_CC <- usedGO(GOData_DD_Up_CC)
allRes_DD_Up_CC <- GenTable(GOData_DD_Up_CC, 
                            classicFisher = resultsFis_Class_DD_Up_CC,
                            weightFisher = resultsFis_Weight_DD_Up_CC,
                            #classicKS = resultsKS_Class_DD_Up_CC,
                            #weightKS = resultsKS_Weight_DD_Up_CC,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_DD_Up_CC))

pVals_Weight_DD_Up_CC <- score(resultsFis_Weight_DD_Up_CC)[score(resultsFis_Weight_DD_Up_CC)<=0.05]

GOData_DD_Up_CC_DE <- termStat(object = GOData_DD_Up_CC, whichGO = names(pVals_Weight_DD_Up_CC))
GOData_DD_Up_CC_DE$DEG <- GOData_DD_Up_CC_DE$Significant
GOData_DD_Up_CC_DE$Term <- Term(rownames(GOData_DD_Up_CC_DE))

ggplot(GOData_DD_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_DD_Up_CC_DE,
            "topGO_weightFS_DD_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## DD Down-regulated

geneList_DD_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_DD_Down$Niben))
names(geneList_DD_Down) <- geneNames_Niben261

GOData_DD_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_DD_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_DD_Down_CC <- runTest(GOData_DD_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_DD_Down_CC <- runTest(GOData_DD_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_DD_Down_CC <- runTest(GOData_DD_Down_CC, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_DD_Down_CC <- runTest(GOData_DD_Down_CC, algorithm = "weight01",
                                        statistic = "fisher")

allGO_DD_Down_CC <- usedGO(GOData_DD_Down_CC)
allRes_DD_Down_CC <- GenTable(GOData_DD_Down_CC, 
                              classicFisher = resultsFis_Class_DD_Down_CC,
                              weightFisher = resultsFis_Weight_DD_Down_CC,
                              #classicKS = resultsKS_Class_DD_Down_CC,
                              #weightKS = resultsKS_Weight_DD_Down_CC,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_DD_Down_CC))

pVals_Weight_DD_Down_CC <- score(resultsFis_Weight_DD_Down_CC)[score(resultsFis_Weight_DD_Down_CC)<=0.05]

GOData_DD_Down_CC_DE <- termStat(object = GOData_DD_Down_CC, whichGO = names(pVals_Weight_DD_Down_CC))
GOData_DD_Down_CC_DE$DEG <- GOData_DD_Down_CC_DE$Significant
GOData_DD_Down_CC_DE$Term <- Term(rownames(GOData_DD_Down_CC_DE))

ggplot(GOData_DD_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_DD_Down_CC_DE,
            "topGO_weightFS_DD_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## EE Up-regulated

geneList_EE_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_EE_Up$Niben))
names(geneList_EE_Up) <- geneNames_Niben261

GOData_EE_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_EE_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_EE_Up_CC <- runTest(GOData_EE_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_EE_Up_CC <- runTest(GOData_EE_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_EE_Up_CC <- runTest(GOData_EE_Up_CC, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_EE_Up_CC <- runTest(GOData_EE_Up_CC, algorithm = "weight01",
                                      statistic = "fisher")

allGO_EE_Up_CC <- usedGO(GOData_EE_Up_CC)
allRes_EE_Up_CC <- GenTable(GOData_EE_Up_CC, 
                            classicFisher = resultsFis_Class_EE_Up_CC,
                            weightFisher = resultsFis_Weight_EE_Up_CC,
                            #classicKS = resultsKS_Class_EE_Up_CC,
                            #weightKS = resultsKS_Weight_EE_Up_CC,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_EE_Up_CC))

pVals_Weight_EE_Up_CC <- score(resultsFis_Weight_EE_Up_CC)[score(resultsFis_Weight_EE_Up_CC)<=0.05]

GOData_EE_Up_CC_DE <- termStat(object = GOData_EE_Up_CC, whichGO = names(pVals_Weight_EE_Up_CC))
GOData_EE_Up_CC_DE$DEG <- GOData_EE_Up_CC_DE$Significant
GOData_EE_Up_CC_DE$Term <- Term(rownames(GOData_EE_Up_CC_DE))

ggplot(GOData_EE_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_EE_Up_CC_DE,
            "topGO_weightFS_EE_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## EE Down-regulated

geneList_EE_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_EE_Down$Niben))
names(geneList_EE_Down) <- geneNames_Niben261

GOData_EE_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_EE_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_EE_Down_CC <- runTest(GOData_EE_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_EE_Down_CC <- runTest(GOData_EE_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_EE_Down_CC <- runTest(GOData_EE_Down_CC, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_EE_Down_CC <- runTest(GOData_EE_Down_CC, algorithm = "weight01",
                                        statistic = "fisher")

allGO_EE_Down_CC <- usedGO(GOData_EE_Down_CC)
allRes_EE_Down_CC <- GenTable(GOData_EE_Down_CC, 
                              classicFisher = resultsFis_Class_EE_Down_CC,
                              weightFisher = resultsFis_Weight_EE_Down_CC,
                              #classicKS = resultsKS_Class_EE_Down_CC,
                              #weightKS = resultsKS_Weight_EE_Down_CC,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_EE_Down_CC))

pVals_Weight_EE_Down_CC <- score(resultsFis_Weight_EE_Down_CC)[score(resultsFis_Weight_EE_Down_CC)<=0.05]

GOData_EE_Down_CC_DE <- termStat(object = GOData_EE_Down_CC, whichGO = names(pVals_Weight_EE_Down_CC))
GOData_EE_Down_CC_DE$DEG <- GOData_EE_Down_CC_DE$Significant
GOData_EE_Down_CC_DE$Term <- Term(rownames(GOData_EE_Down_CC_DE))

ggplot(GOData_EE_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_EE_Down_CC_DE,
            "topGO_weightFS_EE_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## YY Up-regulated

geneList_YY_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_YY_Up$Niben))
names(geneList_YY_Up) <- geneNames_Niben261

GOData_YY_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_YY_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_YY_Up_CC <- runTest(GOData_YY_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_YY_Up_CC <- runTest(GOData_YY_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_YY_Up_CC <- runTest(GOData_YY_Up_CC, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_YY_Up_CC <- runTest(GOData_YY_Up_CC, algorithm = "weight01",
                                      statistic = "fisher")

allGO_YY_Up_CC <- usedGO(GOData_YY_Up_CC)
allRes_YY_Up_CC <- GenTable(GOData_YY_Up_CC, 
                            classicFisher = resultsFis_Class_YY_Up_CC,
                            weightFisher = resultsFis_Weight_YY_Up_CC,
                            #classicKS = resultsKS_Class_YY_Up_CC,
                            #weightKS = resultsKS_Weight_YY_Up_CC,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_YY_Up_CC))

pVals_Weight_YY_Up_CC <- score(resultsFis_Weight_YY_Up_CC)[score(resultsFis_Weight_YY_Up_CC)<=0.05]

GOData_YY_Up_CC_DE <- termStat(object = GOData_YY_Up_CC, whichGO = names(pVals_Weight_YY_Up_CC))
GOData_YY_Up_CC_DE$DEG <- GOData_YY_Up_CC_DE$Significant
GOData_YY_Up_CC_DE$Term <- Term(rownames(GOData_YY_Up_CC_DE))

ggplot(GOData_YY_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_YY_Up_CC_DE,
            "topGO_weightFS_YY_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## YY Down-regulated

geneList_YY_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_YY_Down$Niben))
names(geneList_YY_Down) <- geneNames_Niben261

GOData_YY_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_YY_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_YY_Down_CC <- runTest(GOData_YY_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_YY_Down_CC <- runTest(GOData_YY_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_YY_Down_CC <- runTest(GOData_YY_Down_CC, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_YY_Down_CC <- runTest(GOData_YY_Down_CC, algorithm = "weight01",
                                        statistic = "fisher")

allGO_YY_Down_CC <- usedGO(GOData_YY_Down_CC)
allRes_YY_Down_CC <- GenTable(GOData_YY_Down_CC, 
                              classicFisher = resultsFis_Class_YY_Down_CC,
                              weightFisher = resultsFis_Weight_YY_Down_CC,
                              #classicKS = resultsKS_Class_YY_Down_CC,
                              #weightKS = resultsKS_Weight_YY_Down_CC,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_YY_Down_CC))

pVals_Weight_YY_Down_CC <- score(resultsFis_Weight_YY_Down_CC)[score(resultsFis_Weight_YY_Down_CC)<=0.05]

GOData_YY_Down_CC_DE <- termStat(object = GOData_YY_Down_CC, whichGO = names(pVals_Weight_YY_Down_CC))
GOData_YY_Down_CC_DE$DEG <- GOData_YY_Down_CC_DE$Significant
GOData_YY_Down_CC_DE$Term <- Term(rownames(GOData_YY_Down_CC_DE))

ggplot(GOData_YY_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_YY_Down_CC_DE,
            "topGO_weightFS_YY_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## ZZ Up-regulated

geneList_ZZ_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_ZZ_Up$Niben))
names(geneList_ZZ_Up) <- geneNames_Niben261

GOData_ZZ_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_ZZ_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_ZZ_Up_CC <- runTest(GOData_ZZ_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_ZZ_Up_CC <- runTest(GOData_ZZ_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_ZZ_Up_CC <- runTest(GOData_ZZ_Up_CC, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_ZZ_Up_CC <- runTest(GOData_ZZ_Up_CC, algorithm = "weight01",
                                      statistic = "fisher")

allGO_ZZ_Up_CC <- usedGO(GOData_ZZ_Up_CC)
allRes_ZZ_Up_CC <- GenTable(GOData_ZZ_Up_CC, 
                            classicFisher = resultsFis_Class_ZZ_Up_CC,
                            weightFisher = resultsFis_Weight_ZZ_Up_CC,
                            #classicKS = resultsKS_Class_ZZ_Up_CC,
                            #weightKS = resultsKS_Weight_ZZ_Up_CC,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_ZZ_Up_CC))

pVals_Weight_ZZ_Up_CC <- score(resultsFis_Weight_ZZ_Up_CC)[score(resultsFis_Weight_ZZ_Up_CC)<=0.05]

GOData_ZZ_Up_CC_DE <- termStat(object = GOData_ZZ_Up_CC, whichGO = names(pVals_Weight_ZZ_Up_CC))
GOData_ZZ_Up_CC_DE$DEG <- GOData_ZZ_Up_CC_DE$Significant
GOData_ZZ_Up_CC_DE$Term <- Term(rownames(GOData_ZZ_Up_CC_DE))

ggplot(GOData_ZZ_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_ZZ_Up_CC_DE,
            "topGO_weightFS_ZZ_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## ZZ Down-regulated

geneList_ZZ_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_ZZ_Down$Niben))
names(geneList_ZZ_Down) <- geneNames_Niben261

GOData_ZZ_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_ZZ_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_ZZ_Down_CC <- runTest(GOData_ZZ_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_ZZ_Down_CC <- runTest(GOData_ZZ_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_ZZ_Down_CC <- runTest(GOData_ZZ_Down_CC, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_ZZ_Down_CC <- runTest(GOData_ZZ_Down_CC, algorithm = "weight01",
                                        statistic = "fisher")

allGO_ZZ_Down_CC <- usedGO(GOData_ZZ_Down_CC)
allRes_ZZ_Down_CC <- GenTable(GOData_ZZ_Down_CC, 
                              classicFisher = resultsFis_Class_ZZ_Down_CC,
                              weightFisher = resultsFis_Weight_ZZ_Down_CC,
                              #classicKS = resultsKS_Class_ZZ_Down_CC,
                              #weightKS = resultsKS_Weight_ZZ_Down_CC,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_ZZ_Down_CC))

pVals_Weight_ZZ_Down_CC <- score(resultsFis_Weight_ZZ_Down_CC)[score(resultsFis_Weight_ZZ_Down_CC)<=0.05]

GOData_ZZ_Down_CC_DE <- termStat(object = GOData_ZZ_Down_CC, whichGO = names(pVals_Weight_ZZ_Down_CC))
GOData_ZZ_Down_CC_DE$DEG <- GOData_ZZ_Down_CC_DE$Significant
GOData_ZZ_Down_CC_DE$Term <- Term(rownames(GOData_ZZ_Down_CC_DE))

ggplot(GOData_ZZ_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_ZZ_Down_CC_DE,
            "topGO_weightFS_ZZ_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

################################################################################################################
## UpSetR diagrams of CC in Cornell TimePoints


Cornell_PCA_Up_CC <- list("AAUp"= c(GOData_AA_Up_CC_DE$Term),
                          "BBUp"= c(GOData_BB_Up_CC_DE$Term),
                          "YYUp"= c(GOData_YY_Up_CC_DE$Term),
                          "ZZUp"= c(GOData_ZZ_Up_CC_DE$Term),
                          "CCUp"= c(GOData_CC_Up_CC_DE$Term),
                          "DDUp"= c(GOData_DD_Up_CC_DE$Term),
                          "EEUp"= c(GOData_EE_Up_CC_DE$Term),
                          "56hUp"= c(GOData_56h_Up_CC_DE$Term))

upset(fromList(Cornell_PCA_Up_CC), sets= rev(c("AAUp", "BBUp", "YYUp", "ZZUp", "CCUp", "DDUp", "EEUp", "56hUp")), order.by = "freq", keep.order = TRUE)

Cornell_PCA_Down_CC <- list("AADown"= c(GOData_AA_Down_CC_DE$Term),
                            "BBDown"= c(GOData_BB_Down_CC_DE$Term),
                            "YYDown"= c(GOData_YY_Down_CC_DE$Term),
                            "ZZDown"= c(GOData_ZZ_Down_CC_DE$Term),
                            "CCDown"= c(GOData_CC_Down_CC_DE$Term),
                            "DDDown"= c(GOData_DD_Down_CC_DE$Term),
                            "EEDown"= c(GOData_EE_Down_CC_DE$Term),
                            "56hDown"= c(GOData_56h_Down_CC_DE$Term))

upset(fromList(Cornell_PCA_Down_CC), sets= rev(c("AADown", "BBDown", "YYDown", "ZZDown", "CCDown", "DDDown", "EEDown", "56hDown")), order.by = "freq", keep.order = TRUE)


#####################################################################################
########################### NICOTIANA BRIARDO #######################################
#####################################################################################

## First, we upload the DEGs files (mapped to Niben261 and analyzed with DESeq2)

setwd("E:/RNAseq/Cornell/06_GeneOntology/Cornell")

## We use the files created in the Cornell_DESeq2.R Script

Cornell_6h_Up <- as.data.frame(resSig_6h_2[resSig_6h_2$log2FoldChange>0,])
Cornell_6h_Up$Niben <- paste0(rownames(Cornell_6h_Up),".1")
Cornell_6h_Down <- as.data.frame(resSig_6h_2[resSig_6h_2$log2FoldChange<0,])
Cornell_6h_Down$Niben <- paste0(rownames(Cornell_6h_Down),".1")

Cornell_24h_Up <- as.data.frame(resSig_24h_2[resSig_24h_2$log2FoldChange>0,])
Cornell_24h_Up$Niben <- paste0(rownames(Cornell_24h_Up),".1")
Cornell_24h_Down <- as.data.frame(resSig_24h_2[resSig_24h_2$log2FoldChange<0,])
Cornell_24h_Down$Niben <- paste0(rownames(Cornell_24h_Down),".1")

Cornell_96h_Up <- as.data.frame(resSig_96h_2[resSig_96h_2$log2FoldChange>0,])
Cornell_96h_Up$Niben <- paste0(rownames(Cornell_96h_Up),".1")
Cornell_96h_Down <- as.data.frame(resSig_96h_2[resSig_96h_2$log2FoldChange<0,])
Cornell_96h_Down$Niben <- paste0(rownames(Cornell_96h_Down),".1")


#############################  BIOLOGICAL PROCESS  ###############################################
#################################################################################################
## 6h Up-regulated

geneList_6h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_6h_Up$Niben))
names(geneList_6h_Up) <- geneNames_Niben261

GOData_6h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_6h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_6h_Up_BP <- runTest(GOData_6h_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_6h_Up_BP <- runTest(GOData_6h_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_6h_Up_BP <- runTest(GOData_6h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_6h_Up_BP <- runTest(GOData_6h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_6h_Up_BP <- usedGO(GOData_6h_Up_BP)
allRes_6h_Up_BP <- GenTable(GOData_6h_Up_BP, 
                             classicFisher = resultsFis_Class_6h_Up_BP,
                             weightFisher = resultsFis_Weight_6h_Up_BP,
                             #classicKS = resultsKS_Class_6h_Up_BP,
                             #weightKS = resultsKS_Weight_6h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_6h_Up_BP))

pVals_Weight_6h_Up_BP <- score(resultsFis_Weight_6h_Up_BP)[score(resultsFis_Weight_6h_Up_BP)<=0.05]

GOData_6h_Up_BP_DE <- termStat(object = GOData_6h_Up_BP, whichGO = names(pVals_Weight_6h_Up_BP))
GOData_6h_Up_BP_DE$DEG <- GOData_6h_Up_BP_DE$Significant
GOData_6h_Up_BP_DE$Term <- Term(rownames(GOData_6h_Up_BP_DE))

ggplot(GOData_6h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_6h_Up_BP_DE,
            "topGO_weightFS_6h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 6h Down-regulated

geneList_6h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_6h_Down$Niben))
names(geneList_6h_Down) <- geneNames_Niben261

GOData_6h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_6h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_6h_Down_BP <- runTest(GOData_6h_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_6h_Down_BP <- runTest(GOData_6h_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_6h_Down_BP <- runTest(GOData_6h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_6h_Down_BP <- runTest(GOData_6h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_6h_Down_BP <- usedGO(GOData_6h_Down_BP)
allRes_6h_Down_BP <- GenTable(GOData_6h_Down_BP, 
                               classicFisher = resultsFis_Class_6h_Down_BP,
                               weightFisher = resultsFis_Weight_6h_Down_BP,
                               #classicKS = resultsKS_Class_6h_Down_BP,
                               #weightKS = resultsKS_Weight_6h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_6h_Down_BP))

pVals_Weight_6h_Down_BP <- score(resultsFis_Weight_6h_Down_BP)[score(resultsFis_Weight_6h_Down_BP)<=0.05]

GOData_6h_Down_BP_DE <- termStat(object = GOData_6h_Down_BP, whichGO = names(pVals_Weight_6h_Down_BP))
GOData_6h_Down_BP_DE$DEG <- GOData_6h_Down_BP_DE$Significant
GOData_6h_Down_BP_DE$Term <- Term(rownames(GOData_6h_Down_BP_DE))

ggplot(GOData_6h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_6h_Down_BP_DE,
            "topGO_weightFS_6h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 24h Up-regulated

geneList_24h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_24h_Up$Niben))
names(geneList_24h_Up) <- geneNames_Niben261

GOData_24h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_24h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_24h_Up_BP <- runTest(GOData_24h_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_24h_Up_BP <- runTest(GOData_24h_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_24h_Up_BP <- runTest(GOData_24h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_24h_Up_BP <- runTest(GOData_24h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_24h_Up_BP <- usedGO(GOData_24h_Up_BP)
allRes_24h_Up_BP <- GenTable(GOData_24h_Up_BP, 
                             classicFisher = resultsFis_Class_24h_Up_BP,
                             weightFisher = resultsFis_Weight_24h_Up_BP,
                             #classicKS = resultsKS_Class_24h_Up_BP,
                             #weightKS = resultsKS_Weight_24h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_24h_Up_BP))

pVals_Weight_24h_Up_BP <- score(resultsFis_Weight_24h_Up_BP)[score(resultsFis_Weight_24h_Up_BP)<=0.05]

GOData_24h_Up_BP_DE <- termStat(object = GOData_24h_Up_BP, whichGO = names(pVals_Weight_24h_Up_BP))
GOData_24h_Up_BP_DE$DEG <- GOData_24h_Up_BP_DE$Significant
GOData_24h_Up_BP_DE$Term <- Term(rownames(GOData_24h_Up_BP_DE))

ggplot(GOData_24h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_24h_Up_BP_DE,
            "topGO_weightFS_24h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 24h Down-regulated

geneList_24h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_24h_Down$Niben))
names(geneList_24h_Down) <- geneNames_Niben261

GOData_24h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_24h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_24h_Down_BP <- runTest(GOData_24h_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_24h_Down_BP <- runTest(GOData_24h_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_24h_Down_BP <- runTest(GOData_24h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_24h_Down_BP <- runTest(GOData_24h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_24h_Down_BP <- usedGO(GOData_24h_Down_BP)
allRes_24h_Down_BP <- GenTable(GOData_24h_Down_BP, 
                               classicFisher = resultsFis_Class_24h_Down_BP,
                               weightFisher = resultsFis_Weight_24h_Down_BP,
                               #classicKS = resultsKS_Class_24h_Down_BP,
                               #weightKS = resultsKS_Weight_24h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_24h_Down_BP))

pVals_Weight_24h_Down_BP <- score(resultsFis_Weight_24h_Down_BP)[score(resultsFis_Weight_24h_Down_BP)<=0.05]

GOData_24h_Down_BP_DE <- termStat(object = GOData_24h_Down_BP, whichGO = names(pVals_Weight_24h_Down_BP))
GOData_24h_Down_BP_DE$DEG <- GOData_24h_Down_BP_DE$Significant
GOData_24h_Down_BP_DE$Term <- Term(rownames(GOData_24h_Down_BP_DE))

ggplot(GOData_24h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_24h_Down_BP_DE,
            "topGO_weightFS_24h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 96h Up-regulated

geneList_96h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_96h_Up$Niben))
names(geneList_96h_Up) <- geneNames_Niben261

GOData_96h_Up_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_96h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_96h_Up_BP <- runTest(GOData_96h_Up_BP, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_96h_Up_BP <- runTest(GOData_96h_Up_BP, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_96h_Up_BP <- runTest(GOData_96h_Up_BP, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_96h_Up_BP <- runTest(GOData_96h_Up_BP, algorithm = "weight01",
                                       statistic = "fisher")

allGO_96h_Up_BP <- usedGO(GOData_96h_Up_BP)
allRes_96h_Up_BP <- GenTable(GOData_96h_Up_BP, 
                             classicFisher = resultsFis_Class_96h_Up_BP,
                             weightFisher = resultsFis_Weight_96h_Up_BP,
                             #classicKS = resultsKS_Class_96h_Up_BP,
                             #weightKS = resultsKS_Weight_96h_Up_BP,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_96h_Up_BP))

pVals_Weight_96h_Up_BP <- score(resultsFis_Weight_96h_Up_BP)[score(resultsFis_Weight_96h_Up_BP)<=0.05]

GOData_96h_Up_BP_DE <- termStat(object = GOData_96h_Up_BP, whichGO = names(pVals_Weight_96h_Up_BP))
GOData_96h_Up_BP_DE$DEG <- GOData_96h_Up_BP_DE$Significant
GOData_96h_Up_BP_DE$Term <- Term(rownames(GOData_96h_Up_BP_DE))

ggplot(GOData_96h_Up_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_96h_Up_BP_DE,
            "topGO_weightFS_96h_Up_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 96h Down-regulated

geneList_96h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_96h_Down$Niben))
names(geneList_96h_Down) <- geneNames_Niben261

GOData_96h_Down_BP <- new("topGOdata", ontology= "BP", allGenes= geneList_96h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_96h_Down_BP <- runTest(GOData_96h_Down_BP, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_96h_Down_BP <- runTest(GOData_96h_Down_BP, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_96h_Down_BP <- runTest(GOData_96h_Down_BP, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_96h_Down_BP <- runTest(GOData_96h_Down_BP, algorithm = "weight01",
                                         statistic = "fisher")

allGO_96h_Down_BP <- usedGO(GOData_96h_Down_BP)
allRes_96h_Down_BP <- GenTable(GOData_96h_Down_BP, 
                               classicFisher = resultsFis_Class_96h_Down_BP,
                               weightFisher = resultsFis_Weight_96h_Down_BP,
                               #classicKS = resultsKS_Class_96h_Down_BP,
                               #weightKS = resultsKS_Weight_96h_Down_BP,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_96h_Down_BP))

pVals_Weight_96h_Down_BP <- score(resultsFis_Weight_96h_Down_BP)[score(resultsFis_Weight_96h_Down_BP)<=0.05]

GOData_96h_Down_BP_DE <- termStat(object = GOData_96h_Down_BP, whichGO = names(pVals_Weight_96h_Down_BP))
GOData_96h_Down_BP_DE$DEG <- GOData_96h_Down_BP_DE$Significant
GOData_96h_Down_BP_DE$Term <- Term(rownames(GOData_96h_Down_BP_DE))

ggplot(GOData_96h_Down_BP_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_96h_Down_BP_DE,
            "topGO_weightFS_96h_Down_BP.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)


################################################################################################################
## UpSetR diagrams of BP in Briardo + Count Matrix based groups of Cornell


Cornell_Briardo_Matrix_UP_BP <- list("22hUp"= c(GOData_22h_Up_BP_DE$Term),
                                     "YYUp"= c(GOData_YY_Up_BP_DE$Term),
                                     "ZZUp"= c(GOData_ZZ_Up_BP_DE$Term),
                                     "40hUp"= c(GOData_40h_Up_BP_DE$Term),
                                     "46hUp"= c(GOData_46h_Up_BP_DE$Term),
                                     "56hUp"= c(GOData_56h_Up_BP_DE$Term),
                                     "96hUp"= c(GOData_96h_Up_BP_DE$Term))

upset(fromList(Cornell_Briardo_Matrix_UP_BP), sets= rev(c("22hUp", "YYUp", "ZZUp", "40hUp", "46hUp","56hUp", "96hUp")), order.by = "freq", keep.order = TRUE)

Cornell_Briardo_Matrix_Down_BP <- list("22hDown"= c(GOData_22h_Down_BP_DE$Term),
                                     "YYDown"= c(GOData_YY_Down_BP_DE$Term),
                                     "ZZDown"= c(GOData_ZZ_Down_BP_DE$Term),
                                     "40hDown"= c(GOData_40h_Down_BP_DE$Term),
                                     "46hDown"= c(GOData_46h_Down_BP_DE$Term),
                                     "56hDown"= c(GOData_56h_Down_BP_DE$Term),
                                     "96hDown"= c(GOData_96h_Down_BP_DE$Term))

upset(fromList(Cornell_Briardo_Matrix_Down_BP), sets= rev(c("22hDown", "YYDown", "ZZDown", "40hDown", "46hDown","56hDown", "96hDown")), order.by = "freq", keep.order = TRUE)

#############################  MOLECULAR FUNCTION  ###############################################
#################################################################################################
## 6h Up-regulated

geneList_6h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_6h_Up$Niben))
names(geneList_6h_Up) <- geneNames_Niben261

GOData_6h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_6h_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_6h_Up_MF <- runTest(GOData_6h_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_6h_Up_MF <- runTest(GOData_6h_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_6h_Up_MF <- runTest(GOData_6h_Up_MF, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_6h_Up_MF <- runTest(GOData_6h_Up_MF, algorithm = "weight01",
                                      statistic = "fisher")

allGO_6h_Up_MF <- usedGO(GOData_6h_Up_MF)
allRes_6h_Up_MF <- GenTable(GOData_6h_Up_MF, 
                            classicFisher = resultsFis_Class_6h_Up_MF,
                            weightFisher = resultsFis_Weight_6h_Up_MF,
                            #classicKS = resultsKS_Class_6h_Up_MF,
                            #weightKS = resultsKS_Weight_6h_Up_MF,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_6h_Up_MF))

pVals_Weight_6h_Up_MF <- score(resultsFis_Weight_6h_Up_MF)[score(resultsFis_Weight_6h_Up_MF)<=0.05]

GOData_6h_Up_MF_DE <- termStat(object = GOData_6h_Up_MF, whichGO = names(pVals_Weight_6h_Up_MF))
GOData_6h_Up_MF_DE$DEG <- GOData_6h_Up_MF_DE$Significant
GOData_6h_Up_MF_DE$Term <- Term(rownames(GOData_6h_Up_MF_DE))

ggplot(GOData_6h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_6h_Up_MF_DE,
            "topGO_weightFS_6h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 6h Down-regulated

geneList_6h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_6h_Down$Niben))
names(geneList_6h_Down) <- geneNames_Niben261

GOData_6h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_6h_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_6h_Down_MF <- runTest(GOData_6h_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_6h_Down_MF <- runTest(GOData_6h_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_6h_Down_MF <- runTest(GOData_6h_Down_MF, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_6h_Down_MF <- runTest(GOData_6h_Down_MF, algorithm = "weight01",
                                        statistic = "fisher")

allGO_6h_Down_MF <- usedGO(GOData_6h_Down_MF)
allRes_6h_Down_MF <- GenTable(GOData_6h_Down_MF, 
                              classicFisher = resultsFis_Class_6h_Down_MF,
                              weightFisher = resultsFis_Weight_6h_Down_MF,
                              #classicKS = resultsKS_Class_6h_Down_MF,
                              #weightKS = resultsKS_Weight_6h_Down_MF,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_6h_Down_MF))

pVals_Weight_6h_Down_MF <- score(resultsFis_Weight_6h_Down_MF)[score(resultsFis_Weight_6h_Down_MF)<=0.05]

GOData_6h_Down_MF_DE <- termStat(object = GOData_6h_Down_MF, whichGO = names(pVals_Weight_6h_Down_MF))
GOData_6h_Down_MF_DE$DEG <- GOData_6h_Down_MF_DE$Significant
GOData_6h_Down_MF_DE$Term <- Term(rownames(GOData_6h_Down_MF_DE))

ggplot(GOData_6h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_6h_Down_MF_DE,
            "topGO_weightFS_6h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 24h Up-regulated

geneList_24h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_24h_Up$Niben))
names(geneList_24h_Up) <- geneNames_Niben261

GOData_24h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_24h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_24h_Up_MF <- runTest(GOData_24h_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_24h_Up_MF <- runTest(GOData_24h_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_24h_Up_MF <- runTest(GOData_24h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_24h_Up_MF <- runTest(GOData_24h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_24h_Up_MF <- usedGO(GOData_24h_Up_MF)
allRes_24h_Up_MF <- GenTable(GOData_24h_Up_MF, 
                             classicFisher = resultsFis_Class_24h_Up_MF,
                             weightFisher = resultsFis_Weight_24h_Up_MF,
                             #classicKS = resultsKS_Class_24h_Up_MF,
                             #weightKS = resultsKS_Weight_24h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_24h_Up_MF))

pVals_Weight_24h_Up_MF <- score(resultsFis_Weight_24h_Up_MF)[score(resultsFis_Weight_24h_Up_MF)<=0.05]

GOData_24h_Up_MF_DE <- termStat(object = GOData_24h_Up_MF, whichGO = names(pVals_Weight_24h_Up_MF))
GOData_24h_Up_MF_DE$DEG <- GOData_24h_Up_MF_DE$Significant
GOData_24h_Up_MF_DE$Term <- Term(rownames(GOData_24h_Up_MF_DE))

ggplot(GOData_24h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_24h_Up_MF_DE,
            "topGO_weightFS_24h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 24h Down-regulated

geneList_24h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_24h_Down$Niben))
names(geneList_24h_Down) <- geneNames_Niben261

GOData_24h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_24h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_24h_Down_MF <- runTest(GOData_24h_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_24h_Down_MF <- runTest(GOData_24h_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_24h_Down_MF <- runTest(GOData_24h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_24h_Down_MF <- runTest(GOData_24h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_24h_Down_MF <- usedGO(GOData_24h_Down_MF)
allRes_24h_Down_MF <- GenTable(GOData_24h_Down_MF, 
                               classicFisher = resultsFis_Class_24h_Down_MF,
                               weightFisher = resultsFis_Weight_24h_Down_MF,
                               #classicKS = resultsKS_Class_24h_Down_MF,
                               #weightKS = resultsKS_Weight_24h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_24h_Down_MF))

pVals_Weight_24h_Down_MF <- score(resultsFis_Weight_24h_Down_MF)[score(resultsFis_Weight_24h_Down_MF)<=0.05]

GOData_24h_Down_MF_DE <- termStat(object = GOData_24h_Down_MF, whichGO = names(pVals_Weight_24h_Down_MF))
GOData_24h_Down_MF_DE$DEG <- GOData_24h_Down_MF_DE$Significant
GOData_24h_Down_MF_DE$Term <- Term(rownames(GOData_24h_Down_MF_DE))

ggplot(GOData_24h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_24h_Down_MF_DE,
            "topGO_weightFS_24h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 96h Up-regulated

geneList_96h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_96h_Up$Niben))
names(geneList_96h_Up) <- geneNames_Niben261

GOData_96h_Up_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_96h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_96h_Up_MF <- runTest(GOData_96h_Up_MF, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_96h_Up_MF <- runTest(GOData_96h_Up_MF, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_96h_Up_MF <- runTest(GOData_96h_Up_MF, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_96h_Up_MF <- runTest(GOData_96h_Up_MF, algorithm = "weight01",
                                       statistic = "fisher")

allGO_96h_Up_MF <- usedGO(GOData_96h_Up_MF)
allRes_96h_Up_MF <- GenTable(GOData_96h_Up_MF, 
                             classicFisher = resultsFis_Class_96h_Up_MF,
                             weightFisher = resultsFis_Weight_96h_Up_MF,
                             #classicKS = resultsKS_Class_96h_Up_MF,
                             #weightKS = resultsKS_Weight_96h_Up_MF,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_96h_Up_MF))

pVals_Weight_96h_Up_MF <- score(resultsFis_Weight_96h_Up_MF)[score(resultsFis_Weight_96h_Up_MF)<=0.05]

GOData_96h_Up_MF_DE <- termStat(object = GOData_96h_Up_MF, whichGO = names(pVals_Weight_96h_Up_MF))
GOData_96h_Up_MF_DE$DEG <- GOData_96h_Up_MF_DE$Significant
GOData_96h_Up_MF_DE$Term <- Term(rownames(GOData_96h_Up_MF_DE))

ggplot(GOData_96h_Up_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_96h_Up_MF_DE,
            "topGO_weightFS_96h_Up_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 96h Down-regulated

geneList_96h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_96h_Down$Niben))
names(geneList_96h_Down) <- geneNames_Niben261

GOData_96h_Down_MF <- new("topGOdata", ontology= "MF", allGenes= geneList_96h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_96h_Down_MF <- runTest(GOData_96h_Down_MF, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_96h_Down_MF <- runTest(GOData_96h_Down_MF, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_96h_Down_MF <- runTest(GOData_96h_Down_MF, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_96h_Down_MF <- runTest(GOData_96h_Down_MF, algorithm = "weight01",
                                         statistic = "fisher")

allGO_96h_Down_MF <- usedGO(GOData_96h_Down_MF)
allRes_96h_Down_MF <- GenTable(GOData_96h_Down_MF, 
                               classicFisher = resultsFis_Class_96h_Down_MF,
                               weightFisher = resultsFis_Weight_96h_Down_MF,
                               #classicKS = resultsKS_Class_96h_Down_MF,
                               #weightKS = resultsKS_Weight_96h_Down_MF,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_96h_Down_MF))

pVals_Weight_96h_Down_MF <- score(resultsFis_Weight_96h_Down_MF)[score(resultsFis_Weight_96h_Down_MF)<=0.05]

GOData_96h_Down_MF_DE <- termStat(object = GOData_96h_Down_MF, whichGO = names(pVals_Weight_96h_Down_MF))
GOData_96h_Down_MF_DE$DEG <- GOData_96h_Down_MF_DE$Significant
GOData_96h_Down_MF_DE$Term <- Term(rownames(GOData_96h_Down_MF_DE))

ggplot(GOData_96h_Down_MF_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_96h_Down_MF_DE,
            "topGO_weightFS_96h_Down_MF.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)


################################################################################################################
## UpSetR diagrams of MF in Briardo + Count Matrix based groups of Cornell


Cornell_Briardo_Matrix_UP_MF <- list("22hUp"= c(GOData_22h_Up_MF_DE$Term),
                                     "YYUp"= c(GOData_YY_Up_MF_DE$Term),
                                     "ZZUp"= c(GOData_ZZ_Up_MF_DE$Term),
                                     "40hUp"= c(GOData_40h_Up_MF_DE$Term),
                                     "46hUp"= c(GOData_46h_Up_MF_DE$Term),
                                     "56hUp"= c(GOData_56h_Up_MF_DE$Term),
                                     "96hUp"= c(GOData_96h_Up_MF_DE$Term))

upset(fromList(Cornell_Briardo_Matrix_UP_MF), sets= rev(c("22hUp", "YYUp", "ZZUp", "40hUp", "46hUp","56hUp", "96hUp")), order.by = "freq", keep.order = TRUE)

Cornell_Briardo_Matrix_Down_MF <- list("22hDown"= c(GOData_22h_Down_MF_DE$Term),
                                       "YYDown"= c(GOData_YY_Down_MF_DE$Term),
                                       "ZZDown"= c(GOData_ZZ_Down_MF_DE$Term),
                                       "40hDown"= c(GOData_40h_Down_MF_DE$Term),
                                       "46hDown"= c(GOData_46h_Down_MF_DE$Term),
                                       "56hDown"= c(GOData_56h_Down_MF_DE$Term),
                                       "96hDown"= c(GOData_96h_Down_MF_DE$Term))

upset(fromList(Cornell_Briardo_Matrix_Down_MF), sets= rev(c("22hDown", "YYDown", "ZZDown", "40hDown", "46hDown","56hDown", "96hDown")), order.by = "freq", keep.order = TRUE)

#############################  CELLULAR COMPONENT  ###############################################
#################################################################################################
## 6h Up-regulated

geneList_6h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_6h_Up$Niben))
names(geneList_6h_Up) <- geneNames_Niben261

GOData_6h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_6h_Up, 
                       annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_6h_Up_CC <- runTest(GOData_6h_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_6h_Up_CC <- runTest(GOData_6h_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_6h_Up_CC <- runTest(GOData_6h_Up_CC, algorithm = "classic", 
                                     statistic = "fisher")
resultsFis_Weight_6h_Up_CC <- runTest(GOData_6h_Up_CC, algorithm = "weight01",
                                      statistic = "fisher")

allGO_6h_Up_CC <- usedGO(GOData_6h_Up_CC)
allRes_6h_Up_CC <- GenTable(GOData_6h_Up_CC, 
                            classicFisher = resultsFis_Class_6h_Up_CC,
                            weightFisher = resultsFis_Weight_6h_Up_CC,
                            #classicKS = resultsKS_Class_6h_Up_CC,
                            #weightKS = resultsKS_Weight_6h_Up_CC,
                            orderBy = "weightFisher",
                            topNodes = length(allGO_6h_Up_CC))

pVals_Weight_6h_Up_CC <- score(resultsFis_Weight_6h_Up_CC)[score(resultsFis_Weight_6h_Up_CC)<=0.05]

GOData_6h_Up_CC_DE <- termStat(object = GOData_6h_Up_CC, whichGO = names(pVals_Weight_6h_Up_CC))
GOData_6h_Up_CC_DE$DEG <- GOData_6h_Up_CC_DE$Significant
GOData_6h_Up_CC_DE$Term <- Term(rownames(GOData_6h_Up_CC_DE))

ggplot(GOData_6h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_6h_Up_CC_DE,
            "topGO_weightFS_6h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 6h Down-regulated

geneList_6h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_6h_Down$Niben))
names(geneList_6h_Down) <- geneNames_Niben261

GOData_6h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_6h_Down, 
                         annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_6h_Down_CC <- runTest(GOData_6h_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_6h_Down_CC <- runTest(GOData_6h_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_6h_Down_CC <- runTest(GOData_6h_Down_CC, algorithm = "classic", 
                                       statistic = "fisher")
resultsFis_Weight_6h_Down_CC <- runTest(GOData_6h_Down_CC, algorithm = "weight01",
                                        statistic = "fisher")

allGO_6h_Down_CC <- usedGO(GOData_6h_Down_CC)
allRes_6h_Down_CC <- GenTable(GOData_6h_Down_CC, 
                              classicFisher = resultsFis_Class_6h_Down_CC,
                              weightFisher = resultsFis_Weight_6h_Down_CC,
                              #classicKS = resultsKS_Class_6h_Down_CC,
                              #weightKS = resultsKS_Weight_6h_Down_CC,
                              orderBy = "weightFisher",
                              topNodes = length(allGO_6h_Down_CC))

pVals_Weight_6h_Down_CC <- score(resultsFis_Weight_6h_Down_CC)[score(resultsFis_Weight_6h_Down_CC)<=0.05]

GOData_6h_Down_CC_DE <- termStat(object = GOData_6h_Down_CC, whichGO = names(pVals_Weight_6h_Down_CC))
GOData_6h_Down_CC_DE$DEG <- GOData_6h_Down_CC_DE$Significant
GOData_6h_Down_CC_DE$Term <- Term(rownames(GOData_6h_Down_CC_DE))

ggplot(GOData_6h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_6h_Down_CC_DE,
            "topGO_weightFS_6h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 24h Up-regulated

geneList_24h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_24h_Up$Niben))
names(geneList_24h_Up) <- geneNames_Niben261

GOData_24h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_24h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_24h_Up_CC <- runTest(GOData_24h_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_24h_Up_CC <- runTest(GOData_24h_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_24h_Up_CC <- runTest(GOData_24h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_24h_Up_CC <- runTest(GOData_24h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_24h_Up_CC <- usedGO(GOData_24h_Up_CC)
allRes_24h_Up_CC <- GenTable(GOData_24h_Up_CC, 
                             classicFisher = resultsFis_Class_24h_Up_CC,
                             weightFisher = resultsFis_Weight_24h_Up_CC,
                             #classicKS = resultsKS_Class_24h_Up_CC,
                             #weightKS = resultsKS_Weight_24h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_24h_Up_CC))

pVals_Weight_24h_Up_CC <- score(resultsFis_Weight_24h_Up_CC)[score(resultsFis_Weight_24h_Up_CC)<=0.05]

GOData_24h_Up_CC_DE <- termStat(object = GOData_24h_Up_CC, whichGO = names(pVals_Weight_24h_Up_CC))
GOData_24h_Up_CC_DE$DEG <- GOData_24h_Up_CC_DE$Significant
GOData_24h_Up_CC_DE$Term <- Term(rownames(GOData_24h_Up_CC_DE))

ggplot(GOData_24h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_24h_Up_CC_DE,
            "topGO_weightFS_24h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 24h Down-regulated

geneList_24h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_24h_Down$Niben))
names(geneList_24h_Down) <- geneNames_Niben261

GOData_24h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_24h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_24h_Down_CC <- runTest(GOData_24h_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_24h_Down_CC <- runTest(GOData_24h_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_24h_Down_CC <- runTest(GOData_24h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_24h_Down_CC <- runTest(GOData_24h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_24h_Down_CC <- usedGO(GOData_24h_Down_CC)
allRes_24h_Down_CC <- GenTable(GOData_24h_Down_CC, 
                               classicFisher = resultsFis_Class_24h_Down_CC,
                               weightFisher = resultsFis_Weight_24h_Down_CC,
                               #classicKS = resultsKS_Class_24h_Down_CC,
                               #weightKS = resultsKS_Weight_24h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_24h_Down_CC))

pVals_Weight_24h_Down_CC <- score(resultsFis_Weight_24h_Down_CC)[score(resultsFis_Weight_24h_Down_CC)<=0.05]

GOData_24h_Down_CC_DE <- termStat(object = GOData_24h_Down_CC, whichGO = names(pVals_Weight_24h_Down_CC))
GOData_24h_Down_CC_DE$DEG <- GOData_24h_Down_CC_DE$Significant
GOData_24h_Down_CC_DE$Term <- Term(rownames(GOData_24h_Down_CC_DE))

ggplot(GOData_24h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_24h_Down_CC_DE,
            "topGO_weightFS_24h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 96h Up-regulated

geneList_96h_Up <- factor(as.integer(geneNames_Niben261 %in% Cornell_96h_Up$Niben))
names(geneList_96h_Up) <- geneNames_Niben261

GOData_96h_Up_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_96h_Up, 
                        annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_96h_Up_CC <- runTest(GOData_96h_Up_CC, algorithm = "weight01", 
#                                     statistic = "ks")
#resultsKS_Class_96h_Up_CC <- runTest(GOData_96h_Up_CC, algorithm = "classic", 
#                                    statistic = "ks")
resultsFis_Class_96h_Up_CC <- runTest(GOData_96h_Up_CC, algorithm = "classic", 
                                      statistic = "fisher")
resultsFis_Weight_96h_Up_CC <- runTest(GOData_96h_Up_CC, algorithm = "weight01",
                                       statistic = "fisher")

allGO_96h_Up_CC <- usedGO(GOData_96h_Up_CC)
allRes_96h_Up_CC <- GenTable(GOData_96h_Up_CC, 
                             classicFisher = resultsFis_Class_96h_Up_CC,
                             weightFisher = resultsFis_Weight_96h_Up_CC,
                             #classicKS = resultsKS_Class_96h_Up_CC,
                             #weightKS = resultsKS_Weight_96h_Up_CC,
                             orderBy = "weightFisher",
                             topNodes = length(allGO_96h_Up_CC))

pVals_Weight_96h_Up_CC <- score(resultsFis_Weight_96h_Up_CC)[score(resultsFis_Weight_96h_Up_CC)<=0.05]

GOData_96h_Up_CC_DE <- termStat(object = GOData_96h_Up_CC, whichGO = names(pVals_Weight_96h_Up_CC))
GOData_96h_Up_CC_DE$DEG <- GOData_96h_Up_CC_DE$Significant
GOData_96h_Up_CC_DE$Term <- Term(rownames(GOData_96h_Up_CC_DE))

ggplot(GOData_96h_Up_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_96h_Up_CC_DE,
            "topGO_weightFS_96h_Up_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)

## 96h Down-regulated

geneList_96h_Down <- factor(as.integer(geneNames_Niben261 %in% Cornell_96h_Down$Niben))
names(geneList_96h_Down) <- geneNames_Niben261

GOData_96h_Down_CC <- new("topGOdata", ontology= "CC", allGenes= geneList_96h_Down, 
                          annot= annFUN.gene2GO, gene2GO= geneID2GO_Niben261)
#resultsKS_Weight_96h_Down_CC <- runTest(GOData_96h_Down_CC, algorithm = "weight01", 
#                                      statistic = "ks")
#resultsKS_Class_96h_Down_CC <- runTest(GOData_96h_Down_CC, algorithm = "classic", 
#                                     statistic = "ks")
resultsFis_Class_96h_Down_CC <- runTest(GOData_96h_Down_CC, algorithm = "classic", 
                                        statistic = "fisher")
resultsFis_Weight_96h_Down_CC <- runTest(GOData_96h_Down_CC, algorithm = "weight01",
                                         statistic = "fisher")

allGO_96h_Down_CC <- usedGO(GOData_96h_Down_CC)
allRes_96h_Down_CC <- GenTable(GOData_96h_Down_CC, 
                               classicFisher = resultsFis_Class_96h_Down_CC,
                               weightFisher = resultsFis_Weight_96h_Down_CC,
                               #classicKS = resultsKS_Class_96h_Down_CC,
                               #weightKS = resultsKS_Weight_96h_Down_CC,
                               orderBy = "weightFisher",
                               topNodes = length(allGO_96h_Down_CC))

pVals_Weight_96h_Down_CC <- score(resultsFis_Weight_96h_Down_CC)[score(resultsFis_Weight_96h_Down_CC)<=0.05]

GOData_96h_Down_CC_DE <- termStat(object = GOData_96h_Down_CC, whichGO = names(pVals_Weight_96h_Down_CC))
GOData_96h_Down_CC_DE$DEG <- GOData_96h_Down_CC_DE$Significant
GOData_96h_Down_CC_DE$Term <- Term(rownames(GOData_96h_Down_CC_DE))

ggplot(GOData_96h_Down_CC_DE, aes(x= Significant/Expected, y= Term))+
  geom_point(aes(size=Significant))

write.table(GOData_96h_Down_CC_DE,
            "topGO_weightFS_96h_Down_CC.txt", sep = "\t", dec = ",",
            quote = FALSE, row.names = TRUE)


################################################################################################################
## UpSetR diagrams of CC in Briardo + Count Matrix based groups of Cornell


Cornell_Briardo_Matrix_UP_CC <- list("22hUp"= c(GOData_22h_Up_CC_DE$Term),
                                     "YYUp"= c(GOData_YY_Up_CC_DE$Term),
                                     "ZZUp"= c(GOData_ZZ_Up_CC_DE$Term),
                                     "40hUp"= c(GOData_40h_Up_CC_DE$Term),
                                     "46hUp"= c(GOData_46h_Up_CC_DE$Term),
                                     "56hUp"= c(GOData_56h_Up_CC_DE$Term),
                                     "96hUp"= c(GOData_96h_Up_CC_DE$Term))

upset(fromList(Cornell_Briardo_Matrix_UP_CC), sets= rev(c("22hUp", "YYUp", "ZZUp", "40hUp", "46hUp","56hUp", "96hUp")), order.by = "freq", keep.order = TRUE)

Cornell_Briardo_Matrix_Down_CC <- list("22hDown"= c(GOData_22h_Down_CC_DE$Term),
                                       "YYDown"= c(GOData_YY_Down_CC_DE$Term),
                                       "ZZDown"= c(GOData_ZZ_Down_CC_DE$Term),
                                       "40hDown"= c(GOData_40h_Down_CC_DE$Term),
                                       "46hDown"= c(GOData_46h_Down_CC_DE$Term),
                                       "56hDown"= c(GOData_56h_Down_CC_DE$Term),
                                       "96hDown"= c(GOData_96h_Down_CC_DE$Term))

upset(fromList(Cornell_Briardo_Matrix_Down_CC), sets= rev(c("22hDown", "YYDown", "ZZDown", "40hDown", "46hDown","56hDown", "96hDown")), order.by = "freq", keep.order = TRUE)



#####################################################################################################################
###############################  UPSET DIAGRAMS BENTHAMIANA - TOMATO   ##############################################

## BIOLOGICAL PROCESS

Cornell_Briardo_Tomato_BP_Up <- list("22hUp"= c(GOData_22h_Up_BP_DE$Term),
                                     "YYUp"= c(GOData_YY_Up_BP_DE$Term),
                                     "ZZUp"= c(GOData_ZZ_Up_BP_DE$Term),
                                     "40hUp"= c(GOData_40h_Up_BP_DE$Term),
                                     "46hUp"= c(GOData_46h_Up_BP_DE$Term),
                                     "56hUp"= c(GOData_56h_Up_BP_DE$Term),
                                     "96hUp"= c(GOData_96h_Up_BP_DE$Term),
                                     "PkUp" = c(GOData_Pk_Up_BP_DE$Term),
                                     "LRUp" = c(GOData_LR_Up_BP_DE$Term),
                                     "RRUp" = c(GOData_RR_Up_BP_DE$Term))

PlotA <- upset(fromList(Cornell_Briardo_Tomato_BP_Up), sets = rev(c("22hUp", "YYUp", "ZZUp", "40hUp", "46hUp", "56hUp", "96hUp", "PkUp", "LRUp", "RRUp")), order.by = "freq", keep.order = TRUE)

pdf("upsetPlot_BP_Up_Cornell_Briardo_Tomato.pdf", width=12, height=8, useDingbats=FALSE)
PlotA
dev.off()


Cornell_Briardo_Tomato_BP_Down <- list("22hDown"= c(GOData_22h_Down_BP_DE$Term),
                                     "YYDown"= c(GOData_YY_Down_BP_DE$Term),
                                     "ZZDown"= c(GOData_ZZ_Down_BP_DE$Term),
                                     "40hDown"= c(GOData_40h_Down_BP_DE$Term),
                                     "46hDown"= c(GOData_46h_Down_BP_DE$Term),
                                     "56hDown"= c(GOData_56h_Down_BP_DE$Term),
                                     "96hDown"= c(GOData_96h_Down_BP_DE$Term),
                                     "PkDown" = c(GOData_Pk_Down_BP_DE$Term),
                                     "LRDown" = c(GOData_LR_Down_BP_DE$Term),
                                     "RRDown" = c(GOData_RR_Down_BP_DE$Term))

PlotB <- upset(fromList(Cornell_Briardo_Tomato_BP_Down), sets = rev(c("22hDown", "YYDown", "ZZDown", "40hDown", "46hDown", "56hDown", "96hDown", "PkDown", "LRDown", "RRDown")), order.by = "freq", keep.order = TRUE)

pdf("upsetPlot_BP_Down_Cornell_Briardo_Tomato.pdf", width=12, height=8, useDingbats=FALSE)
PlotB
dev.off()

## MOLECULAR FUNCTION

Cornell_Briardo_Tomato_MF_Up <- list("22hUp"= c(GOData_22h_Up_MF_DE$Term),
                                     "YYUp"= c(GOData_YY_Up_MF_DE$Term),
                                     "ZZUp"= c(GOData_ZZ_Up_MF_DE$Term),
                                     "40hUp"= c(GOData_40h_Up_MF_DE$Term),
                                     "46hUp"= c(GOData_46h_Up_MF_DE$Term),
                                     "56hUp"= c(GOData_56h_Up_MF_DE$Term),
                                     "96hUp"= c(GOData_96h_Up_MF_DE$Term),
                                     "PkUp" = c(GOData_Pk_Up_MF_DE$Term),
                                     "LRUp" = c(GOData_LR_Up_MF_DE$Term),
                                     "RRUp" = c(GOData_RR_Up_MF_DE$Term))

PlotC <- upset(fromList(Cornell_Briardo_Tomato_MF_Up), sets = rev(c("22hUp", "YYUp", "ZZUp", "40hUp", "46hUp", "56hUp", "96hUp", "PkUp", "LRUp", "RRUp")), order.by = "freq", keep.order = TRUE)

pdf("upsetPlot_MF_Up_Cornell_Briardo_Tomato.pdf", width=12, height=8, useDingbats=FALSE)
PlotC
dev.off()


Cornell_Briardo_Tomato_MF_Down <- list("22hDown"= c(GOData_22h_Down_MF_DE$Term),
                                       "YYDown"= c(GOData_YY_Down_MF_DE$Term),
                                       "ZZDown"= c(GOData_ZZ_Down_MF_DE$Term),
                                       "40hDown"= c(GOData_40h_Down_MF_DE$Term),
                                       "46hDown"= c(GOData_46h_Down_MF_DE$Term),
                                       "56hDown"= c(GOData_56h_Down_MF_DE$Term),
                                       "96hDown"= c(GOData_96h_Down_MF_DE$Term),
                                       "PkDown" = c(GOData_Pk_Down_MF_DE$Term),
                                       "LRDown" = c(GOData_LR_Down_MF_DE$Term),
                                       "RRDown" = c(GOData_RR_Down_MF_DE$Term))

PlotD <- upset(fromList(Cornell_Briardo_Tomato_MF_Down), sets = rev(c("22hDown", "YYDown", "ZZDown", "40hDown", "46hDown", "56hDown", "96hDown", "PkDown", "LRDown", "RRDown")), order.by = "freq", keep.order = TRUE)

pdf("upsetPlot_MF_Down_Cornell_Briardo_Tomato.pdf", width=12, height=8, useDingbats=FALSE)
PlotD
dev.off()

## CELLULAR COMPONENT

Cornell_Briardo_Tomato_CC_Up <- list("22hUp"= c(GOData_22h_Up_CC_DE$Term),
                                     "YYUp"= c(GOData_YY_Up_CC_DE$Term),
                                     "ZZUp"= c(GOData_ZZ_Up_CC_DE$Term),
                                     "40hUp"= c(GOData_40h_Up_CC_DE$Term),
                                     "46hUp"= c(GOData_46h_Up_CC_DE$Term),
                                     "56hUp"= c(GOData_56h_Up_CC_DE$Term),
                                     "96hUp"= c(GOData_96h_Up_CC_DE$Term),
                                     "PkUp" = c(GOData_Pk_Up_CC_DE$Term),
                                     "LRUp" = c(GOData_LR_Up_CC_DE$Term),
                                     "RRUp" = c(GOData_RR_Up_CC_DE$Term))

PlotE <- upset(fromList(Cornell_Briardo_Tomato_CC_Up), sets = rev(c("22hUp", "YYUp", "ZZUp", "40hUp", "46hUp", "56hUp", "96hUp", "PkUp", "LRUp", "RRUp")), order.by = "freq", keep.order = TRUE)

pdf("upsetPlot_CC_Up_Cornell_Briardo_Tomato.pdf", width=12, height=8, useDingbats=FALSE)
PlotE
dev.off()


Cornell_Briardo_Tomato_CC_Down <- list("22hDown"= c(GOData_22h_Down_CC_DE$Term),
                                       "YYDown"= c(GOData_YY_Down_CC_DE$Term),
                                       "ZZDown"= c(GOData_ZZ_Down_CC_DE$Term),
                                       "40hDown"= c(GOData_40h_Down_CC_DE$Term),
                                       "46hDown"= c(GOData_46h_Down_CC_DE$Term),
                                       "56hDown"= c(GOData_56h_Down_CC_DE$Term),
                                       "96hDown"= c(GOData_96h_Down_CC_DE$Term),
                                       "PkDown" = c(GOData_Pk_Down_CC_DE$Term),
                                       "LRDown" = c(GOData_LR_Down_CC_DE$Term),
                                       "RRDown" = c(GOData_RR_Down_CC_DE$Term))

PlotF <- upset(fromList(Cornell_Briardo_Tomato_CC_Down), sets = rev(c("22hDown", "YYDown", "ZZDown", "40hDown", "46hDown", "56hDown", "96hDown", "PkDown", "LRDown", "RRDown")), order.by = "freq", keep.order = TRUE)

pdf("upsetPlot_CC_Down_Cornell_Briardo_Tomato.pdf", width=12, height=8, useDingbats=FALSE)
PlotF
dev.off()


#############################################################################################################
## WITHOUT TOMATO


## BIOLOGICAL PROCESS

Cornell_Briardo_Tomato_BP_Up_2 <- list("1st"= c(GOData_22h_Up_BP_DE$Term),
                                     "2nd"= c(GOData_YY_Up_BP_DE$Term),
                                     "3rd"= c(GOData_ZZ_Up_BP_DE$Term),
                                     "4th"= c(GOData_40h_Up_BP_DE$Term),
                                     "5th"= c(GOData_46h_Up_BP_DE$Term),
                                     "6th"= c(GOData_56h_Up_BP_DE$Term),
                                     "7th"= c(GOData_96h_Up_BP_DE$Term))

Plot2 <- upset(fromList(Cornell_Briardo_Tomato_BP_Up_2), sets = rev(c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th")), order.by = "freq", keep.order = TRUE, sets.bar.color=rev(c('#59cac0', '#87c8a8', '#a6c58f', '#bec376', '#d3c05b', '#e4bd3c', '#f4ba00')), shade.color=rev(c('#59cac0', '#a6c58f', '#d3c05b', '#f4ba00')))

pdf("upsetPlot_BP_Up_Cornell_Briardo.pdf", width=12, height=8, useDingbats=FALSE)
Plot2
dev.off()


Cornell_Briardo_Tomato_BP_Down_2 <- list("1st"= c(GOData_22h_Down_BP_DE$Term),
                                     "2nd"= c(GOData_YY_Down_BP_DE$Term),
                                     "3rd"= c(GOData_ZZ_Down_BP_DE$Term),
                                     "4th"= c(GOData_40h_Down_BP_DE$Term),
                                     "5th"= c(GOData_46h_Down_BP_DE$Term),
                                     "6th"= c(GOData_56h_Down_BP_DE$Term),
                                     "7th"= c(GOData_96h_Down_BP_DE$Term))

Plot3 <- upset(fromList(Cornell_Briardo_Tomato_BP_Down_2), sets = rev(c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th")), order.by = "freq", keep.order = TRUE, sets.bar.color=rev(c('#59cac0', '#87c8a8', '#a6c58f', '#bec376', '#d3c05b', '#e4bd3c', '#f4ba00')), shade.color=rev(c('#59cac0', '#a6c58f', '#d3c05b', '#f4ba00')))

pdf("upsetPlot_BP_Down_Cornell_Briardo.pdf", width=12, height=8, useDingbats=FALSE)
Plot3
dev.off()


## MOLECULAR FUNCTION

Cornell_Briardo_Tomato_MF_Up_2 <- list("1st"= c(GOData_22h_Up_MF_DE$Term),
                                     "2nd"= c(GOData_YY_Up_MF_DE$Term),
                                     "3rd"= c(GOData_ZZ_Up_MF_DE$Term),
                                     "4th"= c(GOData_40h_Up_MF_DE$Term),
                                     "5th"= c(GOData_46h_Up_MF_DE$Term),
                                     "6th"= c(GOData_56h_Up_MF_DE$Term),
                                     "7th"= c(GOData_96h_Up_MF_DE$Term))

Plot4 <- upset(fromList(Cornell_Briardo_Tomato_MF_Up_2), sets = rev(c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th")), order.by = "freq", keep.order = TRUE, sets.bar.color=rev(c('#59cac0', '#87c8a8', '#a6c58f', '#bec376', '#d3c05b', '#e4bd3c', '#f4ba00')), shade.color=rev(c('#59cac0', '#a6c58f', '#d3c05b', '#f4ba00')))

pdf("upsetPlot_MF_Up_Cornell_Briardo.pdf", width=12, height=8, useDingbats=FALSE)
Plot4
dev.off()


Cornell_Briardo_Tomato_MF_Down_2 <- list("1st"= c(GOData_22h_Down_MF_DE$Term),
                                     "2nd"= c(GOData_YY_Down_MF_DE$Term),
                                     "3rd"= c(GOData_ZZ_Down_MF_DE$Term),
                                     "4th"= c(GOData_40h_Down_MF_DE$Term),
                                     "5th"= c(GOData_46h_Down_MF_DE$Term),
                                     "6th"= c(GOData_56h_Down_MF_DE$Term),
                                     "7th"= c(GOData_96h_Down_MF_DE$Term))

Plot5 <- upset(fromList(Cornell_Briardo_Tomato_MF_Down_2), sets = rev(c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th")), order.by = "freq", keep.order = TRUE, sets.bar.color=rev(c('#59cac0', '#87c8a8', '#a6c58f', '#bec376', '#d3c05b', '#e4bd3c', '#f4ba00')), shade.color=rev(c('#59cac0', '#a6c58f', '#d3c05b', '#f4ba00')))

pdf("upsetPlot_MF_Down_Cornell_Briardo.pdf", width=12, height=8, useDingbats=FALSE)
Plot5
dev.off()

## CELLULAR COMPONENT

Cornell_Briardo_Tomato_CC_Up_2 <- list("1st"= c(GOData_22h_Up_CC_DE$Term),
                                     "2nd"= c(GOData_YY_Up_CC_DE$Term),
                                     "3rd"= c(GOData_ZZ_Up_CC_DE$Term),
                                     "4th"= c(GOData_40h_Up_CC_DE$Term),
                                     "5th"= c(GOData_46h_Up_CC_DE$Term),
                                     "6th"= c(GOData_56h_Up_CC_DE$Term),
                                     "7th"= c(GOData_96h_Up_CC_DE$Term))

Plot6 <- upset(fromList(Cornell_Briardo_Tomato_CC_Up_2), sets = rev(c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th")), order.by = "freq", keep.order = TRUE, sets.bar.color=rev(c('#59cac0', '#87c8a8', '#a6c58f', '#bec376', '#d3c05b', '#e4bd3c', '#f4ba00')), shade.color=rev(c('#59cac0', '#a6c58f', '#d3c05b', '#f4ba00')))

pdf("upsetPlot_CC_Up_Cornell_Briardo.pdf", width=12, height=8, useDingbats=FALSE)
Plot6
dev.off()


Cornell_Briardo_Tomato_CC_Down_2 <- list("1st"= c(GOData_22h_Down_CC_DE$Term),
                                     "2nd"= c(GOData_YY_Down_CC_DE$Term),
                                     "3rd"= c(GOData_ZZ_Down_CC_DE$Term),
                                     "4th"= c(GOData_40h_Down_CC_DE$Term),
                                     "5th"= c(GOData_46h_Down_CC_DE$Term),
                                     "6th"= c(GOData_56h_Down_CC_DE$Term),
                                     "7th"= c(GOData_96h_Down_CC_DE$Term))

Plot7 <- upset(fromList(Cornell_Briardo_Tomato_CC_Down_2), sets = rev(c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th")), order.by = "freq", keep.order = TRUE, sets.bar.color=rev(c('#59cac0', '#87c8a8', '#a6c58f', '#bec376', '#d3c05b', '#e4bd3c', '#f4ba00')), shade.color=rev(c('#59cac0', '#a6c58f', '#d3c05b', '#f4ba00')))

pdf("upsetPlot_CC_Down_Cornell_Briardo.pdf", width=12, height=8, useDingbats=FALSE)
Plot7
dev.off()
