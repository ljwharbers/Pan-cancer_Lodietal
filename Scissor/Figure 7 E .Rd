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
################################################################################################


df <- readRDS("GitHub.4_DataPerSample.rds")
head( df)

            tNEG tpos BULK    Module Patient totalscissor tumortype              Outcome  fraction   Subtype
BC_BK11...1   56   96  HCC POSModule BC_BK11          152        BC Positive correlation 0.6315789 POSModule
BC_BK12...2   52  100  HCC POSModule BC_BK12          152        BC Positive correlation 0.6578947 POSModule
BC_BK13...3    3    1  HCC POSModule BC_BK13            4        BC Positive correlation 0.2500000 POSModule
BC_BK14...4    7   12  HCC POSModule BC_BK14           19        BC Positive correlation 0.6315789 POSModule
BC_BK15...5  122  193  HCC POSModule BC_BK15          315        BC Positive correlation 0.6126984 POSModule
BC_BK16...6   10   53  HCC POSModule BC_BK16           63        BC Positive correlation 0.8412698 POSModule

##Select sampletypes in module

TLS10 <- c( "CD4+ TFH", "Memory GC-indep","Breg", "Plasmablast", "IFN Mac", "IgG mature", "IgG immature", "IgA mature","mQuiescDC", "IgA immature")
POSModule <- c("CD4+ TREG", "mRegDC", "Mono-like Mac", "Neutrophils", "LAM2","Inflam Mac" , "AXL_DC", "CD4+ TH1", "CD8+ TEX", "Proliferating" , "Lymphatic")



DF <- df [df$Subtype %in% TLS10, ]
L <- split(DF,DF$BULK)
ll <- list()
for (i in 1:length(L)) 
{
dfs <- L[[i]]
Module <- cbind(
tapply(dfs$tNEG,dfs$Patient, sum),
tapply(dfs$tpos,dfs$Patient, sum))
Module <- as.data.frame(as.matrix(Module))
colnames(Module) <- c("tNEG","tpos")
Module$BULK <- unique(dfs$BULK)
ll[[i]] <- Module
}

TLS10 <- do.call("rbind", ll)
TLS10$Module <- "TLS10"


table(POSModule  %in%df$Subtype)

DF <- df [df$Subtype %in% POSModule , ]

## Get BULK from SET ...
L <- split(DF,DF$BULK)
#Should we keep this seperate per BULK dataset. 
ll <- list()
for (i in 1:length(L)) 
{
dfs <- L[[i]]
Module <- cbind(
tapply(dfs$tNEG,dfs$Patient, sum),
tapply(dfs$tpos,dfs$Patient, sum))
Module <- as.data.frame(as.matrix(Module))
colnames(Module) <- c("tNEG","tpos")
Module$BULK <- unique(dfs$BULK)
ll[[i]] <- Module
}

POSModule <- do.call("rbind", ll)
POSModule$Module <- "POSModule"


f<- rbind(POSModule, TLS10)
f$Patient <- rownames(f)
f$totalscissor <- f$tNEG + f$tpos  
fp <- f
fp$Outcome <- "Positive correlation"
fn <- f
fn$Outcome <- "Negative correlation"
df <- bind_rows(fp,fn)


df$fraction <- "zero"
df$fraction[which(df$Outcome == "Negative correlation")]  <- df$tNEG[which(df$Outcome == "Negative correlation")] / df$totalscissor[which(df$Outcome == "Negative correlation")]
df$fraction[which(df$Outcome == "Positive correlation")]  <- df$tpos[which(df$Outcome == "Positive correlation")] / df$totalscissor[which(df$Outcome == "Positive correlation")]
df$fraction <- as.numeric(df$fraction)
df$Outcome <-factor(df$Outcome, levels=c("Negative correlation", "Positive correlation" ))



##########################
##########################
################# Figure 7 E
##########################

df$Subtype <- df$Module
MAJOR <- "Modules"
options(ggplot2.discrete.colour = c("midnightblue","red2"))
options(repr.plot.width=21, repr.plot.height=7)
set.seed(123)

# Selecting samples >9 cells
FILTER=9
D <- df[which (df$totalscissor > FILTER), ] 

p2 <- ggplot(D , aes(x=Subtype, y=fraction, color=Outcome)) +
geom_boxplot(outlier.shape=NA)+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_compare_means(method = "t.test", paired = F, label.y=1, label="p.signif") + ggtitle (paste0("Only subtypes with >",FILTER," scissor cells")) +  geom_point(data = D[which (D$BULK=="HCC" ),], aes(shape = "HCC"), position = position_dodge(width=0.9)) +
   geom_point(data = D[which (D$BULK=="NSCLC" ),], aes(shape = "NSCLC"), position = position_dodge(width=0.8)) +
     geom_point(data = D[which (D$BULK=="OAK" ),], aes(shape = "OAK"), position = position_dodge(width=0.7))  +
	   geom_point(data = D[which (D$BULK=="mUC"),], aes(shape = "mUC"), position = position_dodge(width=0.6)) +   
	     guides(fill = guide_legend(override.aes = list(shape = c( NA, NA))))
p2

	 library(gridExtra)
	 pdf(paste0("BoxplotsPerPatient.all.",FILTER,".",MAJOR,".pdf"), width=10, height=6)
p2
	 dev.off()



############################################
################ Supplementary Figure 7 E
############################################
######## Per PanCancer tumortype

f$tumortype <- "NA"
f$tumortype [grep ("^HCC", rownames(f))] <- "HCC"
f$tumortype [grep ("^HNSCC", rownames(f))] <- "HNSCC"
f$tumortype [grep ("^OV", rownames(f))] <- "OV"
f$tumortype [grep ("^MEL", rownames(f))] <- "MEL"
f$tumortype [grep ("^NSCLC", rownames(f))] <- "NSCLC"
f$tumortype [grep ("^GBM", rownames(f))] <- "GBM"
f$tumortype [grep ("^CC", rownames(f))] <- "CC"
f$tumortype [grep ("^BC", rownames(f))] <- "BC"
f$tumortype [grep ("^CRC", rownames(f))] <- "CRC"
f$tumortype [grep ("^CC3", rownames(f))] <- "endometrial"

fp <- f
fp$Outcome <- "Positive correlation"
fn <- f
fn$Outcome <- "Negative correlation"
df <- bind_rows(fp,fn)

df$fraction <- "zero"
df$fraction[which(df$Outcome == "Negative correlation")]  <- df$tNEG[which(df$Outcome == "Negative correlation")] / df$totalscissor[which(df$Outcome == "Negative correlation")]
df$fraction[which(df$Outcome == "Positive correlation")]  <- df$tpos[which(df$Outcome == "Positive correlation")] / df$totalscissor[which(df$Outcome == "Positive correlation")]
df$fraction <- as.numeric(df$fraction)
df$Outcome <-factor(df$Outcome, levels=c("Negative correlation", "Positive correlation" ))


df$Subtype <- df$Module
MAJOR <- "Modules"
options(ggplot2.discrete.colour = c("midnightblue","red2"))
options(repr.plot.width=21, repr.plot.height=7)


# Selecting samples >9 cells
FILTER=9
D <- df[which (df$totalscissor > FILTER), ] 
D$tumortype <- gsub("endometrial", "CC", D$tumortype)

p1 <- ggplot(D , aes(x=Subtype, y=fraction, color=Outcome)) +
geom_boxplot(outlier.shape=NA)+
geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.3), size = 0.4, alpha = 0.7) + facet_wrap(~tumortype) +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) ) + stat_compare_means(method = "t.test", paired = F, label.y=1, label="p.signif") + ggtitle (paste0("Only subtypes with >",FILTER," scissor cells, based on 4 bulk datasets"))

pdf(paste0("BoxplotsPerPatientperTumorTypeModule.pdf"), width=21, height=14)
print(p1)
dev.off()







