
print.noquote("bboeckx, V2, check results carefully")
print.noquote("load also Library (dplyr)") 
EnsembleToGenes <- function(object) {
test <- read.csv("/data/refgenomes/Genes10X/Final.Ensemble.csv", row.names=1,stringsAsFactors=F)
UNI <- read.csv("/data/refgenomes/Genes10X/uniqueV3.0_GEX.csv",row.names=1,stringsAsFactors=F)
object2 <-  object [(UNI[,1]),]


rownames(object2) <- plyr::mapvalues(x = rownames(object2), from = test$V1, to = test$V2)
  return(object2)
}

