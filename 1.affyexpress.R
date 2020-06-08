rm(list = ls())
options(stringsAsFactors = F)

#download raw data
library(ArrayExpress)
AE = getAE("E-MTAB-5563",path="./RAW files/")
# save(AE,file = "AE.Rdata")
rawset= ae2bioc(mageFiles = AE)

#RMA
library(oligo)
AEsetnorm = oligo::rma(rawset)
exprSet = exprs(AEsetnorm)
boxplot(exprSet,outline = F) #check normalizatin
pdata.all = pData(AEsetnorm)
pdata = pdata.all[,c(3,6:9)]
colnames(exprSet) = unlist(lapply(colnames(exprSet), function(x) {
  strsplit(x,"\\.")[[1]][1]
}))
rownames(pdata) = unlist(lapply(rownames(pdata), function(x) {
  strsplit(x,"\\.")[[1]][1]
}))
boxplot(exprSet)
pdata = pdata[colnames(exprSet),]
pdata$Characteristics.disease. = c(rep("typeI",5),rep("typeI_control",6),
                                   rep("typeII_control",6),rep("typeII",6),
                                   rep("insulres_control",6),rep("insulres",6))

colnames(pdata) = c("strain","genotype","phenotype","diet","disease")
pdata$type = c(rep("typeI",5),rep("control",6),
               rep("control",6),rep("typeII",6),
               rep("control",6),rep("insulres",6))

save(pdata,pdata.all,exprSet,file = "eset.Rdata")
