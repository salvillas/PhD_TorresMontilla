library(tidyverse)
library(UpSetR)
library(grid)
library(topGO)
library(VennDiagram)
library(RColorBrewer)
library(ggvenn)
library(ggplot2)
library(ggpubr)

for (file in grep("weightFS", value = TRUE, x= grep("_2.txt", list.files("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/"), value = TRUE))){
  nam <- read.table(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/",file), sep = "\t", header = FALSE)
  name <- str_sub(file, 1, -7)
  terms <- as.data.frame(Term(nam$V1))
  terms <- terms[!is.na(terms)]
  write.table(terms, paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/", name, "_3.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

topGO_weightFS_BP_Up <- grep("weightFS", value = TRUE, x = grep("Up_BP_2.txt", list.files("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/"), value = TRUE))
topGO_weightFS_BP_Down <- grep("weightFS", value = TRUE, x = grep("Down_BP_2.txt", list.files("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/"), value = TRUE))
topGO_weightFS_MF_Up <- grep("weightFS", value = TRUE, x = grep("Up_MF_2.txt", list.files("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/"), value = TRUE))
topGO_weightFS_MF_Down <- grep("weightFS", value = TRUE, x = grep("Down_MF_2.txt", list.files("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/"), value = TRUE))
topGO_weightFS_CC_Up <- grep("weightFS", value = TRUE, x = grep("Up_CC_2.txt", list.files("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/"), value = TRUE))
topGO_weightFS_CC_Down <- grep("weightFS", value = TRUE, x = grep("Down_CC_2.txt", list.files("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/"), value = TRUE))

Cornell_TimePoints <- c("22h", "25h", "28h", "34h", "37h", "40h", "46h", "56h", "96h")

Upset_plot_GO <- function(Category){
  lista <- list()
  for (file in grep(paste(Cornell_TimePoints, collapse = "|"),get(paste0("topGO_weightFS_", Category)), value = TRUE)){
    nam <- read.table(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/", file), 
                      sep = "\t", header = FALSE)
    name <- str_sub(file, 16, 18)
    lista[[name]] <- nam[,1]
  }
  return(upset(fromList(lista), sets = rev(names(lista)), order.by = "freq", keep.order = TRUE))
}

Upset_list_GO <- function(Category){
  lista <- list()
  for (file in grep(paste(Cornell_TimePoints, collapse = "|"),get(paste0("topGO_weightFS_", Category)), value = TRUE)){
    nam <- read.table(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/", file), 
                      sep = "\t", header = TRUE)
    name <- str_sub(file, 16, -5)
    lista[[name]] <- nam[,1]
  }
  return(lista)
}


DEGs_Cornell <- grep("filtered_2.txt", x = list.files("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana"), value = TRUE)


Upset_plot_DEG <- function(Category= c("Up", "Down")){
  lista <- list()
  for (file in grep(paste(Cornell_TimePoints, collapse = "|"),DEGs_Cornell, value = TRUE)){
    nam <- read.table(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/", file), 
                      sep = "\t", header = TRUE)
    if (Category=="Up"){
      nam <- nam[nam[,2]>0,]
    }
    else {
      nam <- nam[nam[,2]<0,]
    }
    name <- str_sub(file, 6, 8)
    lista[[name]] <- nam[,1]
  }
  return(upset(fromList(lista), sets = rev(names(lista)), order.by = "freq", keep.order = TRUE))
}

DEGs_Up_plot <- Upset_plot_DEG("Up")
pdf(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/", "upsetPlot_DEGs_Up_TimePoints_Cornell_96h.pdf"),
    height = 6, width = 10)
DEGs_Up_plot; grid.text("Up",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

DEGs_Down_plot <- Upset_plot_DEG("Down")
pdf(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/", "UpsetPlot_DEGs_Down_TimePoints_Cornell_96h.pdf"),
    height = 6, width = 10)
DEGs_Down_plot; grid.text("Down",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

BP_Up_plot <- Upset_plot_GO("BP_Up")
pdf(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/UpSetR Diagrams/", "upsetPlot_BP_Up_TimePoints_Cornell.pdf"),
    height = 6, width = 10)
BP_Up_plot; grid.text("BP Up",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()
BP_Down_plot <- Upset_plot_GO("BP_Down")
pdf(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/UpSetR Diagrams/", "upsetPlot_BP_Down_TimePoints_Cornell.pdf"),
    height = 6, width = 10)
BP_Down_plot; grid.text("BP Down",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

MF_Up_plot <- Upset_plot_GO("MF_Up")
pdf(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/UpSetR Diagrams/", "upsetPlot_MF_Up_TimePoints_Cornell.pdf"),
    height = 6, width = 10)
MF_Up_plot; grid.text("MF Up",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()
MF_Down_plot <- Upset_plot_GO("MF_Down")
pdf(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/UpSetR Diagrams/", "upsetPlot_MF_Down_TimePoints_Cornell.pdf"),
    height = 6, width = 10)
MF_Down_plot; grid.text("MF Down",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

CC_Up_plot <- Upset_plot_GO("CC_Up")
pdf(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/UpSetR Diagrams/", "upsetPlot_CC_Up_TimePoints_Cornell.pdf"),
    height = 6, width = 10)
CC_Up_plot; grid.text("CC Up",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()
CC_Down_plot <- Upset_plot_GO("CC_Down")
pdf(paste0("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/UpSetR Diagrams/", "upsetPlot_CC_Down_TimePoints_Cornell.pdf"),
    height = 6, width = 10)
CC_Down_plot; grid.text("CC Down",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

BP_Up_list <- Upset_list_GO("BP_Up")
BP_Down_list <- Upset_list_GO("BP_Down")
MF_Up_list <- Upset_list_GO("MF_Up")
MF_Down_list <- Upset_list_GO("MF_Down")
CC_Up_list <- Upset_list_GO("CC_Up")
CC_Down_list <- Upset_list_GO("CC_Down")


## Comparison 56h vs 96h

myCol <- brewer.pal(4, "Pastel2")

Venn_diagram <- function (Category){
  Category_split <- str_split(Category, pattern = "_")
  list <- get(paste(Category, "list", sep = "_"))[c(paste("56h", Category_split[[1]][2], Category_split[[1]][1], sep = "_"),
                                                    paste("96h", Category_split[[1]][2], Category_split[[1]][1], sep = "_"))]
  venn.diagram(list, 
               filename = paste0("Venn_56h_96h_", Category, ".png"),
               output= TRUE,
               category.names = c(paste("56h", Category_split[[1]][1], Category_split[[1]][2], sep = " "), 
                                  paste("96h", Category_split[[1]][1], Category_split[[1]][2], sep = " ")),
               imagetype = "png",
               height = 480,
               width = 480,
               resolution = 300,
               compression = "lzw",
               lwd= 2,
               lty= 'blank',
               fill= c("#B3E2CD", "#FDCDAC"),
               cex= 0.6,
               fontface= "bold",
               fontfamily= "sans",
               cat.cex = 0.6,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.pos = c(-20, 20),
               cat.dist = c(0.055, 0.055),
               cat.fontfamily = "sans")
}
Venn_diagram("BP_Up")
Venn_diagram("BP_Down")
Venn_diagram("MF_Up")
Venn_diagram("MF_Down")


DEGs_56h_Up <- read.table("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/pcrtB56h_vs_GFP56h_filtered_2.txt",
                          sep = "\t", header = TRUE, dec = ",") %>% dplyr::filter(pcrtB56h_vs_GFP56h_filtered > 0) %>% dplyr::select(1)
DEGs_56h_Down <- read.table("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/pcrtB56h_vs_GFP56h_filtered_2.txt",
                          sep = "\t", header = TRUE, dec = ",") %>% dplyr::filter(pcrtB56h_vs_GFP56h_filtered < 0) %>% dplyr::select(1)
DEGs_96h_Up <- read.table("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/pcrtB96h_vs_GFP96h_filtered_2.txt",
                          sep = "\t", header = TRUE, dec = ",") %>% dplyr::filter(pcrtB96h_vs_GFP96h_filtered > 0) %>% dplyr::select(1)
DEGs_96h_Down <- read.table("/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/pcrtB96h_vs_GFP96h_filtered_2.txt",
                            sep = "\t", header = TRUE, dec = ",") %>% dplyr::filter(pcrtB96h_vs_GFP96h_filtered < 0) %>% dplyr::select(1)
DEGs_56h_vs_96h <- list("56h_Up" = DEGs_56h_Up$Name,
                        "56h_Down" = DEGs_56h_Down$Name,
                        "96h_Up" = DEGs_96h_Up$Name,
                        "96h_Down" = DEGs_96h_Down$Name)
venn.diagram(DEGs_56h_vs_96h, 
             filename = "/Users/salva/GoogleDrive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/UpSetR Diagrams/Venn_56h_96h_DEGs.png",
             output= TRUE,
             category.names = c("56h Up", "56h Down", "96h Up", "96h Down"),
             imagetype = "png",
             height = 620,
             width = 620,
             resolution = 300,
             compression = "lzw",
             lwd= 2,
             lty= 'blank',
             fill= myCol,
             cex= 0.6,
             fontface= "bold",
             fontfamily= "sans",
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             #cat.pos = c(-20, 20),
             #cat.dist = c(0.055, 0.055),
             cat.fontfamily = "sans")


########################## This is a way I found on Internet to obtain a list with the elements from different groups in the Venn
Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

Get_elements_Upset <- function(list){
  combs <- unlist(lapply(1:length(list), 
                         function(j) combn(names(list), j, simplify = FALSE)),
                  recursive = FALSE)
  names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
  
  elements <- lapply(combs, function(i) Setdiff(list[i], list[setdiff(names(list), i)]))
  return(elements)
}

BP_Up_list_elements <- Get_elements_Upset(BP_Up_list)



###### NEW VENN DIAGRAMS THESIS

BP_56h_Up <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_56h_Up_BP_2.txt",
                        header = FALSE, sep = "\t")
BP_56h_Down <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_56h_Down_BP_2.txt",
                        header = FALSE, sep = "\t")
BP_96h_Up <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_96h_Up_BP_2.txt",
                        header = FALSE, sep = "\t")
BP_96h_Down <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_96h_Down_BP_2.txt",
                          header = FALSE, sep = "\t")

BP_Venn <- list("56h Up"= BP_56h_Up$V1,
                  "56h Down"= BP_56h_Down$V1,
                  "96h Up"= BP_96h_Up$V1,
                  "96h Down"= BP_96h_Down$V1)

BP <- ggvenn(BP_Venn,
             fill_color = c("#9b8872","#3a2c19", "#fba882","#d95f02"),
             stroke_size = 0.5, set_name_size = 4)

MF_56h_Up <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_56h_Up_MF_2.txt",
                        header = FALSE, sep = "\t")
MF_56h_Down <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_56h_Down_MF_2.txt",
                          header = FALSE, sep = "\t")
MF_96h_Up <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_96h_Up_MF_2.txt",
                        header = FALSE, sep = "\t")
MF_96h_Down <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_96h_Down_MF_2.txt",
                          header = FALSE, sep = "\t")


MF_Venn <- list("56h Up"= MF_56h_Up$V1,
                "56h Down"= MF_56h_Down$V1,
                "96h Up"= MF_96h_Up$V1,
                "96h Down"= MF_96h_Down$V1)

MF <- ggvenn(MF_Venn,
             fill_color = c("#9b8872","#3a2c19", "#fba882","#d95f02"),
             stroke_size = 0.5, set_name_size = 4)

CC_56h_Up <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_56h_Up_CC_2.txt",
                        header = FALSE, sep = "\t")
CC_56h_Down <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_56h_Down_CC_2.txt",
                          header = FALSE, sep = "\t")
CC_96h_Up <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_96h_Up_CC_2.txt",
                        header = FALSE, sep = "\t")
CC_96h_Down <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/06_GeneOntology/Cornell/2nd/topGO_weightFS_96h_Down_CC_2.txt",
                          header = FALSE, sep = "\t")


CC_Venn <- list("56h Up"= CC_56h_Up$V1,
                "56h Down"= CC_56h_Down$V1,
                "96h Up"= CC_96h_Up$V1,
                "96h Down"= CC_96h_Down$V1)

CC <- ggvenn(CC_Venn,
             fill_color = c("#9b8872","#3a2c19", "#fba882","#d95f02"),
             stroke_size = 0.5, set_name_size = 4)


DEG_56h <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/pcrtB56h_vs_GFP56h_filtered_2.txt",
                        header = TRUE, sep = "\t")
DEG_56h_Up <- DEG_56h[DEG_56h[,2]>0,]
DEG_56h_Down <- DEG_56h[DEG_56h[,2]<0,]

DEG_96h <- read.table("/Users/salva/Google Drive/My Drive/CRAG/RNAseq/Milan/RNAseq_Milan/DEGs_Transcript/Nicotiana/pcrtB96h_vs_GFP96h_filtered_2.txt",
                      header = TRUE, sep = "\t")
DEG_96h_Up <- DEG_96h[DEG_96h[,2]>0,]
DEG_96h_Down <- DEG_96h[DEG_96h[,2]<0,]

DEG_Venn <- list("56h Up"= DEG_56h_Up$Name,
                "56h Down"= DEG_56h_Down$Name,
                "96h Up"= DEG_96h_Up$Name,
                "96h Down"= DEG_96h_Down$Name)

DEG <- ggvenn(DEG_Venn,
             fill_color = c("#9b8872","#3a2c19", "#fba882","#d95f02"),
             stroke_size = 0.5, set_name_size = 4)

ggarrange(DEG, BP, MF, CC, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))


### Groups of Manuel interest (events (i, ii, iv); 22h, 25-28h, 46-56h)

Event1_Up <- as.data.frame(Term(c("GO:0009734",
                                  "GO:0000160",
                                  "GO:0009435",
                                  "GO:0015689",
                                  "GO:0043622",
                                  "GO:0007010",
                                  "GO:0006811",
                                  "GO:1902476",
                                  "GO:0006270",
                                  "GO:0052865")))
Event1_Up$GOterm <- rownames(Event1_Up)
colnames(Event1_Up) <- c("GOterm","Event1_Up")
write.table(Event1_Up, "/Users/salva/Google Drive/My Drive/CRAG/Tesis/Chapter 2/Tables/BP_UpSetR_Event1_Up.txt",
          row.names = FALSE, sep = "\t", quote = FALSE)

Event2_Up <- as.data.frame(Term(c("GO:0048255",
                                 "GO:0070899",
                                 "GO:0006562",
                                 "GO:0006048",
                                 "GO:0070143",
                                 "GO:1903826",
                                 "GO:0046294",
                                 "GO:0015813",
                                 "GO:0052746",
                                 "GO:0071158",
                                 "GO:0010019",
                                 "GO:1901962",
                                 "GO:0006567",
                                 "GO:0032774",
                                 "GO:0010044",
                                 "GO:0032515",
                                 "GO:1901703",
                                 "GO:0010847",
                                 "GO:1905786",
                                 "GO:0009663",
                                 "GO:0015808",
                                 "GO:0015812",
                                 "GO:0006633",
                                 "GO:0015949",
                                 "GO:1903401",
                                 "GO:0042545",
                                 "GO:0009432")))
Event2_Up$GOterm <- rownames(Event2_Up)
colnames(Event2_Up) <- c("GOterm","Event2_Up")
write.table(Event2_Up, "/Users/salva/Google Drive/My Drive/CRAG/Tesis/Chapter 2/Tables/BP_UpSetR_Event2_Up.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

Event4_Up <- as.data.frame(Term(c("GO:0080142",
                                  "GO:0051762",
                                  "GO:0042256",
                                  "GO:0034204",
                                  "GO:0019441",
                                  "GO:0030001",
                                  "GO:0046513",
                                  "GO:0046777",
                                  "GO:0006952",
                                  "GO:0007178",
                                  "GO:0002229",
                                  "GO:0006598",
                                  "GO:0019722",
                                  "GO:0030187",
                                  "GO:0071456",
                                  "GO:0015914",
                                  "GO:0048544",
                                  "GO:0080156")))
Event4_Up$GOterm <- rownames(Event4_Up)
colnames(Event4_Up) <- c("GOterm","Event4_Up")
write.table(Event4_Up, "/Users/salva/Google Drive/My Drive/CRAG/Tesis/Chapter 2/Tables/BP_UpSetR_Event4_Up.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

Event1_Down <- as.data.frame(Term(c("GO:0019430",
                                    "GO:0080142",
                                    "GO:0009694",
                                    "GO:0048268",
                                    "GO:0048278",
                                    "GO:0007186",
                                    "GO:0046513",
                                    "GO:0006527",
                                    "GO:0046470",
                                    "GO:0009058",
                                    "GO:0006355",
                                    "GO:0006096",
                                    "GO:0007188",
                                    "GO:0050790",
                                    "GO:0022900",
                                    "GO:0031204",
                                    "GO:0006102",
                                    "GO:0006694",
                                    "GO:0032502")))
Event1_Down$GOterm <- rownames(Event1_Down)
colnames(Event1_Down) <- c("GOterm","Event1_Down")
write.table(Event1_Down, "/Users/salva/Google Drive/My Drive/CRAG/Tesis/Chapter 2/Tables/BP_UpSetR_Event1_Down.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

Event2_Down <- as.data.frame(Term(c("GO:0046777",
                                    "GO:0006515",
                                    "GO:0010142",
                                    "GO:0042908",
                                    "GO:0042742",
                                    "GO:0010952",
                                    "GO:0006075",
                                    "GO:0005992",
                                    "GO:2000072",
                                    "GO:0010777",
                                    "GO:0006435",
                                    "GO:0010155",
                                    "GO:0009626",
                                    "GO:0035494",
                                    "GO:0015867",
                                    "GO:0019419")))
Event2_Down$GOterm <- rownames(Event2_Down)
colnames(Event2_Down) <- c("GOterm","Event2_Down")
write.table(Event2_Down, "/Users/salva/Google Drive/My Drive/CRAG/Tesis/Chapter 2/Tables/BP_UpSetR_Event2_Down.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)


Event4_Down <- as.data.frame(Term(c("GO:0010236",
                                    "GO:0008610",
                                    "GO:0005978",
                                    "GO:0006817",
                                    "GO:0098506",
                                    "GO:0019252",
                                    "GO:0048229",
                                    "GO:0034219",
                                    "GO:0046854",
                                    "GO:0006021",
                                    "GO:0006857",
                                    "GO:0009641",
                                    "GO:0009640")))
Event4_Down$GOterm <- rownames(Event4_Down)
colnames(Event4_Down) <- c("GOterm","Event4_Down")
write.table(Event4_Down, "/Users/salva/Google Drive/My Drive/CRAG/Tesis/Chapter 2/Tables/BP_UpSetR_Event4_Down.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)


