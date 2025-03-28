library(ggplot2)
library(tidyverse)
library(Seurat)
library(plyr)
library(metafor)
library(gridExtra)
library(scales)
library(dplyr)
library(stringr)
library(readxl)
library(reshape)
library(ggpubr)


set.seed(123)


################################################################################################
################################################################################################
################################################################################################

DF <- readRDS(file="4_DataPerPhenotype.rds")

> head(DF)
         Subtype tPOS tNEG  cPOS  cNEG                   SET TYPE                 TYPE_SET     SubtypeS TotallCells tPOS_NEG tNS
441 Naive mature 1118 1391 10417 10004 Scissor_OAK_COX_00001   BC BC_Scissor_OAK_COX_00001 Naive mature        2511     2509   2
443  Memory IgM+ 2705 3134 10417 10004 Scissor_OAK_COX_00001   BC BC_Scissor_OAK_COX_00001  Memory IgM+        5849     5839  10
445         Breg 1092  567 10417 10004 Scissor_OAK_COX_00001   BC BC_Scissor_OAK_COX_00001         Breg        1659     1659   0
449   IgG mature 3002 2123 10417 10004 Scissor_OAK_COX_00001   BC BC_Scissor_OAK_COX_00001   IgG mature        5131     5125   6
451 IgA immature  283  373 10417 10004 Scissor_OAK_COX_00001   BC BC_Scissor_OAK_COX_00001 IgA immature         656      656   0
453  Memory IgM-  212  238 10417 10004 Scissor_OAK_COX_00001   BC BC_Scissor_OAK_COX_00001  Memory IgM-         451      450   1
    TotallCellsMajor TotallCells_No_tPOS TotallCells_No_tNEG specific_cPOS specific_TotallCellsMajor BULK
441            20445                1393                1120          9299                     17934  OAK
443            20445                3144                2715          7712                     14596  OAK
445            20445                 567                1092          9325                     18786  OAK
449            20445                2129                3008          7415                     15314  OAK
451            20445                 373                 283         10134                     19789  OAK
453            20445                 239                 213         10205                     19994  OAK
> 

TLS10 <- c( "CD4+ TFH", "Memory GC-indep","Breg", "Plasmablast", "IFN Mac", "IgG mature", "IgG immature", "IgA mature","mQuiescDC", "IgA immature")
POSModule <- c("CD4+ TREG", "mRegDC", "Mono-like Mac", "Neutrophils", "LAM2","Inflam Mac" , "AXL_DC", "CD4+ TH1", "CD8+ TEX", "Proliferating" , "Lymphatic") # We don't have results for ECs so no lymphatic.


############################################################
## Add the Immune hubs 
## Run two times: one time for each hub


SELECT <- "TLS10"
select <- TLS10

SELECT <- "POSModule"
select <- POSModule

O <- DF
dm1 <- O [which (O$SubtypeS %in% select), ]
dm1b <- unique(O [which (O$SubtypeS %in% select),c("cPOS", "cNEG","BULK") ])
cPOS <- as.data.frame(as.matrix(tapply(dm1b$cPOS, dm1b$BULK, sum)))
cNEG <- as.data.frame(as.matrix(tapply(dm1b$cNEG, dm1b$BULK, sum)))
tPOS <- as.data.frame(as.matrix(tapply(dm1$tPOS, dm1$BULK,sum)))
tNEG <- as.data.frame(as.matrix(tapply(dm1$tNEG, dm1$BULK,sum)))
TotallCells <- as.data.frame(as.matrix(tapply(dm1$TotallCells, dm1$BULK,sum)))
dm <- as.data.frame(cbind(tPOS,tNEG,TotallCells, cPOS, cNEG))


colnames(dm) <- c("tPOS", "tNEG","TotallCells","cPOS", "cNEG")
dm$BULK <- rownames( dm)
dm$Subtype <- SELECT
dm$SubtypeS <- SELECT
dm$SET <- rownames(dm)
dm$colorV2b<- colM

dm$tPOS_NEG <- dm$tPOS+ dm$tNEG
dm$tNS <- dm$TotallCells - dm$tPOS_NEG
dm$TotallCellsMajor <- NA

dm$TotallCells_No_tPOS <- dm$TotallCells -dm$tPOS
dm$TotallCells_No_tNEG <- dm$TotallCells -dm$tNEG
dm$specific_cPOS <-dm$cPOS  -dm$tPOS
dm$specific_TotallCellsMajor <- NA

colnames(O) %in% colnames(dm)
L <- levels(O$SubtypeS)
L <- c(L,SELECT)

mod <- bind_rows(O,dm)
mod$SubtypeS <- factor(mod$SubtypeS, levels=L)
O <- mod
DF <- O


############################################################
#### Calculate enrichment

dat1 <- escalc(measure="OR", ai=tPOS, bi=tNEG, ci=cPOS, di=cNEG, data=DF, digits=4)
dat1<- summary.escalc(dat1,digits=5)
head(dat1)


################################################################################################
################################################################################################
## plot Enrichment Figure 7D; Supplementary Figure 7D
################################################################################################
# We calculate the average
DATA <-	dat1 
	N <- data.frame(matrix(ncol = length( colnames(DATA)), nrow = length(unique(DATA$SubtypeS))))
	colnames(DATA) -> colnames(N)
	N$yi <- as.data.frame(as.matrix(tapply(DATA$yi, DATA$SubtypeS,mean ))) [,1]
	N$SubtypeS <- names(tapply(DATA$yi, DATA$SubtypeS,mean ))
	N$BULK <- "Average"
dat1 <- 	rbind(DATA,N)


## Figures 7 D , Supplemenatry Figure 7 D
NAME <- "Average"
	print.noquote(NAME)
DatO <- dat1[which (dat1$BULK==NAME),]

p1<-ggplot(data=DatO, aes(x=reorder(SubtypeS, yi), y=yi, fill=yi)) +
  geom_bar(stat="identity")+
  coord_flip() +
    ggtitle (paste0("Enrichment analysis ", NAME))+
	scale_fill_gradient2(low='dodgerblue4', mid='lavenderblush', high='red4', midpoint=0, , limits = c(-2.5, 2.5))

p2 <- p1 +  geom_point(data = dat1[which (dat1$BULK=="HCC" ),], aes(shape = "HCC"), position = position_dodge(width = 0.9)) +
   geom_point(data = dat1[which (dat1$BULK=="NSCLC" ),], aes(shape = "NSCLC"), position = position_dodge(width = 0.9)) +
     geom_point(data = dat1[which (dat1$BULK=="OAK" ),], aes(shape = "OAK"), position = position_dodge(width = 0.9))  +
	   geom_point(data = dat1[which (dat1$BULK=="mUC"),], aes(shape = "mUC"), position = position_dodge(width = 0.9))
	   
################################################################################################
################################################################################################
################################################################################################
################################################################################################


