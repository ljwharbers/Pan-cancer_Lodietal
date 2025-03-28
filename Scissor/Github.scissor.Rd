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


Bulkdata <- readRDS("Bulkdata.rds")
Bulkmetadata <- readRDS("Bulkmetadata.rds")
Object <- readRDS("TCell.rds") # SeuratObject


#The scissor analysis was applied to the major immune celltypes (TC, BC, MM and DC) with 4 bulk datasets (mUC, HCC, NSCLC(poplar), NSCLC(OAK)). 

Bulkdata[1:5,1:5]
         EA.177c304419e7 EA.4ebd09823c3e EA.a395dd19348f EA.0c979151be41
A1BG          0.71726033       0.7739183     3.817665219       1.9847290
A1BG-AS1      0.99114380       1.1666580     6.851115000       2.0569430
A1CF          0.01856896       0.0000000     0.009103179       0.0000000
A2M         304.69370000     692.5841000    94.960700000    1165.6540000
A2M-AS1       0.44128530       0.5681256     0.283938400       0.8013332


head(Bulkmetadata)
                     time status
EA.177c304419e7 5.5196205      1
EA.4ebd09823c3e 1.5113247      1
EA.a395dd19348f 0.6899526      1
EA.0c979151be41 4.2054252      1
EA.a2714bb09f23 1.4456149      1
EA.1ac6c6502028 2.7269554      1


##alpha values used

BC mUC 0.001
DC mUC 0.01
MM mUC 0.0625
TC mUC 0.001
BC NSCLC 0.001
DC NSCLC 0.0625
MM NSCLC 0.0625
TC NSCLC 0.01
BC OAK 0.0001
DC OAK 0.0001
MM OAK 0.001
TC OAK 0.0625
BC HCC 0.001
DC HCC 0.001
MM HCC 0.0001
TC HCC 0.0001



#The analysis below is run for every 

res_cox<- Scissor(bulk_dataset=as.matrix(Bulkdata), 
                  sc_dataset=Object, 
                  Bulkmetadata, 
                  alpha = alpha,  
                  family = "cox") 
          print.noquote("done")
       print.noquote(Sys.time())
saveRDS(res_cox, file= paste0("Scissor_Result_COX.rds"))

# Results for each separate analysis are merged in the metadata.
meta <- Object@meta.data
meta$Results_Scissor_HCC <- "not_significant"
meta$Results_Scissor_HCC [which(rownames(meta) %in% res_cox$Scissor_pos)]  <- "Scissor_pos"
meta$Results_Scissor_HCC [which(rownames(meta) %in% res_cox$Scissor_neg)]  <- "Scissor_neg"


