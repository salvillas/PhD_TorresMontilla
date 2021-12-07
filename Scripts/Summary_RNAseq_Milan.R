library(tidyverse)
library(DESeq2)
library(ggplot2)
library(gplots)
library(genefilter)
library(GenomicRanges)
library(plyr)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(Biobase)
library(gplots)
library(graphics)
library(grDevices)
library(survival)

## First we load the file with TPM from all genes

TPM_All <- read.table("TPM_all_samples.txt", 
                      sep = " ", header = TRUE, na.strings = "")

## Then we read the annotation and we include it to the previous file

Niben261_annotation <- read_delim("Niben261_annotation.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

TPM_All_Annotated <- merge(Niben261_annotation, TPM_All, by="Name", all = TRUE)

## I load the 6h file filtered and the file with All FC in this comparison

FC_Filtered_6h <- read.table("DEGs_Transcript/pcrtB6h_vs_GFP6h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_6h_Annotated <- merge(FC_Filtered_6h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_6h_Annotated, "DEGs_Transcript/pcrtB6h_vs_GFP6h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

FC_All_6h <- read.table("DEGs_Transcript/pcrtB6h_vs_GFP6h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_6h_Annotated <- merge(FC_All_6h, Niben261_annotation, by="Name", all.x = TRUE)


## I load the 22h file filtered and the file with All FC in this comparison

FC_Filtered_22h <- read.table("DEGs_Transcript/pcrtB22h_vs_GFP22h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_22h_Annotated <- merge(FC_Filtered_22h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_22h_Annotated, "DEGs_Transcript/pcrtB22h_vs_GFP22h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(FC_Filtered_6h, FC_Filtered_22h, by="Name", all = TRUE)


FC_All_22h <- read.table("DEGs_Transcript/pcrtB22h_vs_GFP22h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_22h_Annotated <- merge(FC_All_22h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(FC_All_6h, FC_All_22h, by="Name", all = TRUE)


## I load the 24h file filtered and the file with All FC in this comparison

FC_Filtered_24h <- read.table("DEGs_Transcript/pcrtB24h_vs_GFP24h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_24h_Annotated <- merge(FC_Filtered_24h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_24h_Annotated, "DEGs_Transcript/pcrtB24h_vs_GFP24h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_24h, by="Name", all = TRUE)

FC_All_24h <- read.table("DEGs_Transcript/pcrtB24h_vs_GFP24h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_24h_Annotated <- merge(FC_All_24h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_24h, by="Name", all = TRUE)


## I load the 25h file filtered and the file with All FC in this comparison

FC_Filtered_25h <- read.table("DEGs_Transcript/pcrtB25h_vs_GFP25h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_25h_Annotated <- merge(FC_Filtered_25h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_25h_Annotated, "DEGs_Transcript/pcrtB25h_vs_GFP25h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_25h, by="Name", all = TRUE)

FC_All_25h <- read.table("DEGs_Transcript/pcrtB25h_vs_GFP25h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_25h_Annotated <- merge(FC_All_25h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_25h, by="Name", all = TRUE)


## I load the AA file filtered and the file with All FC in this comparison

FC_Filtered_AA <- read.table("DEGs_Transcript/pcrtBAA_vs_GFPAA_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_AA_Annotated <- merge(FC_Filtered_AA, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_AA_Annotated, "DEGs_Transcript/pcrtBAA_vs_GFPAA_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_AA, by="Name", all = TRUE)

FC_All_AA <- read.table("DEGs_Transcript/pcrtBAA_vs_GFPAA_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_AA_Annotated <- merge(FC_All_AA, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_AA, by="Name", all = TRUE)



## I load the BB file filtered and the file with All FC in this comparison

FC_Filtered_BB <- read.table("DEGs_Transcript/pcrtBBB_vs_GFPBB_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_BB_Annotated <- merge(FC_Filtered_BB, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_BB_Annotated, "DEGs_Transcript/pcrtBBB_vs_GFPBB_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_BB, by="Name", all = TRUE)

FC_All_BB <- read.table("DEGs_Transcript/pcrtBBB_vs_GFPBB_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_BB_Annotated <- merge(FC_All_BB, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_BB, by="Name", all = TRUE)


## I load the YY file filtered and the file with All FC in this comparison

FC_Filtered_YY <- read.table("DEGs_Transcript/pcrtBYY_vs_GFPYY_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_YY_Annotated <- merge(FC_Filtered_YY, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_YY_Annotated, "DEGs_Transcript/pcrtBYY_vs_GFPYY_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_YY, by="Name", all = TRUE)

FC_All_YY <- read.table("DEGs_Transcript/pcrtBYY_vs_GFPYY_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_YY_Annotated <- merge(FC_All_YY, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_YY, by="Name", all = TRUE)


## I load the 28h file filtered and the file with All FC in this comparison

FC_Filtered_28h <- read.table("DEGs_Transcript/pcrtB28h_vs_GFP28h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_28h_Annotated <- merge(FC_Filtered_28h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_28h_Annotated, "DEGs_Transcript/pcrtB28h_vs_GFP28h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_28h, by="Name", all = TRUE)

FC_All_28h <- read.table("DEGs_Transcript/pcrtB28h_vs_GFP28h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_28h_Annotated <- merge(FC_All_28h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_28h, by="Name", all = TRUE)


## I load the 34h file filtered and the file with All FC in this comparison

FC_Filtered_34h <- read.table("DEGs_Transcript/pcrtB34h_vs_GFP34h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_34h_Annotated <- merge(FC_Filtered_34h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_34h_Annotated, "DEGs_Transcript/pcrtB34h_vs_GFP34h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_34h, by="Name", all = TRUE)

FC_All_34h <- read.table("DEGs_Transcript/pcrtB34h_vs_GFP34h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_34h_Annotated <- merge(FC_All_34h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_34h, by="Name", all = TRUE)


## I load the 37h file filtered and the file with All FC in this comparison

FC_Filtered_37h <- read.table("DEGs_Transcript/pcrtB37h_vs_GFP37h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_37h_Annotated <- merge(FC_Filtered_37h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_37h_Annotated, "DEGs_Transcript/pcrtB37h_vs_GFP37h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_37h, by="Name", all = TRUE)

FC_All_37h <- read.table("DEGs_Transcript/pcrtB37h_vs_GFP37h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_37h_Annotated <- merge(FC_All_37h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_37h, by="Name", all = TRUE)


## I load the ZZ file filtered and the file with All FC in this comparison

FC_Filtered_ZZ <- read.table("DEGs_Transcript/pcrtBZZ_vs_GFPZZ_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_ZZ_Annotated <- merge(FC_Filtered_ZZ, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_ZZ_Annotated, "DEGs_Transcript/pcrtBZZ_vs_GFPZZ_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_ZZ, by="Name", all = TRUE)

FC_All_ZZ <- read.table("DEGs_Transcript/pcrtBZZ_vs_GFPZZ_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_ZZ_Annotated <- merge(FC_All_ZZ, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_ZZ, by="Name", all = TRUE)


## I load the CC file filtered and the file with All FC in this comparison

FC_Filtered_CC <- read.table("DEGs_Transcript/pcrtBCC_vs_GFPCC_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_CC_Annotated <- merge(FC_Filtered_CC, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_CC_Annotated, "DEGs_Transcript/pcrtBCC_vs_GFPCC_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_CC, by="Name", all = TRUE)

FC_All_CC <- read.table("DEGs_Transcript/pcrtBCC_vs_GFPCC_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_CC_Annotated <- merge(FC_All_CC, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_CC, by="Name", all = TRUE)


## I load the 40h file filtered and the file with All FC in this comparison

FC_Filtered_40h <- read.table("DEGs_Transcript/pcrtB40h_vs_GFP40h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_40h_Annotated <- merge(FC_Filtered_40h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_40h_Annotated, "DEGs_Transcript/pcrtB40h_vs_GFP40h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_40h, by="Name", all = TRUE)

FC_All_40h <- read.table("DEGs_Transcript/pcrtB40h_vs_GFP40h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_40h_Annotated <- merge(FC_All_40h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_40h, by="Name", all = TRUE)


## I load the DD file filtered and the file with All FC in this comparison

FC_Filtered_DD <- read.table("DEGs_Transcript/pcrtBDD_vs_GFPDD_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_DD_Annotated <- merge(FC_Filtered_DD, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_DD_Annotated, "DEGs_Transcript/pcrtBDD_vs_GFPDD_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_DD, by="Name", all = TRUE)

FC_All_DD <- read.table("DEGs_Transcript/pcrtBDD_vs_GFPDD_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_DD_Annotated <- merge(FC_All_DD, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_DD, by="Name", all = TRUE)


## I load the 46h file filtered and the file with All FC in this comparison

FC_Filtered_46h <- read.table("DEGs_Transcript/pcrtB46h_vs_GFP46h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_46h_Annotated <- merge(FC_Filtered_46h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_46h_Annotated, "DEGs_Transcript/pcrtB46h_vs_GFP46h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_46h, by="Name", all = TRUE)

FC_All_46h <- read.table("DEGs_Transcript/pcrtB46h_vs_GFP46h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_46h_Annotated <- merge(FC_All_46h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_46h, by="Name", all = TRUE)


## I load the 56h file filtered and the file with All FC in this comparison

FC_Filtered_56h <- read.table("DEGs_Transcript/pcrtB56h_vs_GFP56h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_56h_Annotated <- merge(FC_Filtered_56h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_56h_Annotated, "DEGs_Transcript/pcrtB56h_vs_GFP56h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_56h, by="Name", all = TRUE)

FC_All_56h <- read.table("DEGs_Transcript/pcrtB56h_vs_GFP56h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_56h_Annotated <- merge(FC_All_56h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_56h, by="Name", all = TRUE)


## I load the EE file filtered and the file with All FC in this comparison

FC_Filtered_EE <- read.table("DEGs_Transcript/pcrtBEE_vs_GFPEE_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_EE_Annotated <- merge(FC_Filtered_EE, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_EE_Annotated, "DEGs_Transcript/pcrtBEE_vs_GFPEE_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_EE, by="Name", all = TRUE)

FC_All_EE <- read.table("DEGs_Transcript/pcrtBEE_vs_GFPEE_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_EE_Annotated <- merge(FC_All_EE, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_EE, by="Name", all = TRUE)


## I load the 96h file filtered and the file with All FC in this comparison

FC_Filtered_96h <- read.table("DEGs_Transcript/pcrtB96h_vs_GFP96h_filtered_2.txt",
                              sep= "\t", header = TRUE, dec=",")
FC_Filtered_96h_Annotated <- merge(FC_Filtered_96h, Niben261_annotation, by="Name", all.x = TRUE)
write.table(FC_Filtered_96h_Annotated, "DEGs_Transcript/pcrtB96h_vs_GFP96h_filtered_Annotated.txt",
            quote = FALSE, dec = ",", row.names = FALSE)

Summary_FC_Filtered <- merge(Summary_FC_Filtered, FC_Filtered_96h, by="Name", all = TRUE)

FC_All_96h <- read.table("DEGs_Transcript/pcrtB96h_vs_GFP96h_All_2.txt",
                         sep = "\t", header = TRUE, dec=",")
FC_All_96h_Annotated <- merge(FC_All_96h, Niben261_annotation, by="Name", all.x = TRUE)

Summary_FC_All <- merge(Summary_FC_All, FC_All_96h, by="Name", all = TRUE)


Summary_FC_Filtered_Annotation <- merge(Summary_FC_Filtered, Niben261_annotation, by="Name", all.x=TRUE)
Summary_FC_All_Annotation <- merge(Summary_FC_All, Niben261_annotation, by="Name", all.x = TRUE)

TPM_All_FC <- merge(TPM_All, Summary_FC_All, by="Name", all = TRUE)
TPM_All_FC <- merge(TPM_All_FC, Summary_FC_Filtered, by="Name", all = TRUE)

TPM_All_FC_Annotated <- merge(Niben261_annotation, TPM_All_FC, by="Name", all = TRUE)


write.table(TPM_All_FC_Annotated, "Summary_AllSamples_TPM_FC.txt",
            quote = FALSE, dec = ",", row.names = FALSE, sep = "\t")


