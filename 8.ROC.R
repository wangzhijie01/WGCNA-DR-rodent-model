rm(list=ls())
options(stringsAsFactors = F)
library(ROCR)
load("WGCNA.Rdata")

head(datTraits)
design=model.matrix(~0+ datTraits$type)
design=data.frame(design,design[,2]+design[,3])
colnames(design)=c(levels(factor(datTraits$type)),"both_type")

blue.hub <- read.csv ("./blue-top30-MCC.csv", header=TRUE, sep=",",skip=1)
blue.hub.gene=blue.hub$Name
blue.hub.gene

magenta.hub <- read.csv ("./magenta-top30-MCC.csv", header=TRUE, sep=",",skip=1)
magenta.hub.gene=magenta.hub$Name
magenta.hub.gene

train.data=data.frame(datExpr,design)


ROC_plot=function(gene="Cxcl6",data=train.data,y="typeII"){
  col=match(gene,colnames(datExpr))
  pred1 <- prediction(data[,col], data[,y])
  perf1 <- performance(pred1,"tpr","fpr")
  auc1 <- performance(pred1,"auc")
  auc1=auc1@y.values[[1]]
  plot(perf1,col="red",main=gene)
  abline(a=0,b=1,col="blue")
  text(0, 1, labels=paste("AUC =", format.pval(auc1,digits = 3)), adj=c(0, .5))
  return(auc1)
}
ROC_plot(gene="Krtap4-3")
ROC_plot(gene="Ifngr2")
  

#blue module
pdf("blue-hub-ROC.pdf",width = 10,height = 15)
par(mfrow = c(6,5))
blue.filter=unlist(lapply(blue.hub.gene,ROC_plot))
names(blue.filter)=blue.hub.gene
blue.filter=data.frame(gene=blue.hub.gene,trainROC=blue.filter)
dev.off()
write.table(blue.filter,"blue-hub-ROC.txt",row.names = F,col.names = F,quote=F)
blue.hub.filter=blue.hub.gene[blue.filter$trainROC>0.8]

#magenta module
pdf("magenta-hub-ROC.pdf",width = 10,height = 15)
par(mfrow = c(6,5))
magenta.filter=unlist(lapply(magenta.hub.gene,ROC_plot,data=train.data,y="typeII"))
magenta.filter=data.frame(gene=magenta.hub.gene,trainROC=magenta.filter)
dev.off()
write.table(magenta.filter,"magenta-hub-ROC.txt",row.names = F,col.names = F,quote=F)
magenta.hub.filter=magenta.hub.gene[magenta.filter$trainROC>0.8]

##validation
if (F) {
  library(GEOquery)
  load("allGSE.Rdata")
  GSE111465.expr=allGSE[[1]]
  GSE28831.expr=allGSE[[2]]
  GSE55389.expr=allGSE[[3]]
  
  GSE111465.pdata=data.frame(typeI=c(rep(0,6),rep(1,6)))
  GSE28831.pdata=pData(GSE28831)[,43]
  GSE28831.pdata=data.frame(typeI=ifelse(GSE28831.pdata=="STZ",1,0))
  GSE55389.pdata=data.frame(typeII=c(rep(1,4),rep(0,4)))
  
  GSE111465.data=data.frame(t(GSE111465.expr),GSE111465.pdata)
  GSE28831.data=data.frame(t(GSE28831.expr),GSE28831.pdata)
  GSE55389.data=data.frame(t(GSE55389.expr),GSE55389.pdata)
  
  
  ROC_plot2=function(gene="Cxcl6",data=train.data,y="typeII"){
    col=match(gene,colnames(data))
    pred1 <- prediction(data[,col], data[,y])
    perf1 <- performance(pred1,"tpr","fpr")
    auc1 <- performance(pred1,"auc")
    auc1=auc1@y.values[[1]]
    plot(perf1,col="red",main=gene)
    abline(a=0,b=1,col="blue")
    text(0, 1, labels=paste("AUC =", format.pval(auc1,digits = 3)), adj=c(0, .5))
    return(auc1)
  }
  
  # GSE55389
  #blue
  GSE55389.blue.filter=blue.hub.filter[blue.hub.filter%in%colnames(GSE55389.data)]
  length(GSE55389.blue.filter)
  pdf("GSE55389-blue-hub-ROC.pdf",width = 10,height = 15)
  par(mfrow = c(6,5))
  GSE55389.blue.ROC=unlist(lapply(GSE55389.blue.filter,ROC_plot2,data=GSE55389.data,y="typeII"))
  dev.off()
  GSE55389.blue.ROC=data.frame(gene=GSE55389.blue.filter,GSE55389.ROC=GSE55389.blue.ROC)
  blue.ROC=merge(blue.filter,GSE55389.blue.ROC,by = "gene",all.x = T)
  write.csv(blue.ROC,"blue.ROC.csv")
  #magenta
  GSE55389.magenta.filter=magenta.hub.filter[magenta.hub.filter%in%colnames(GSE55389.data)]
  length(GSE55389.magenta.filter)
  pdf("GSE55389-magenta-hub-ROC.pdf",width = 10,height = 15)
  par(mfrow = c(6,5))
  GSE55389.magenta.ROC=unlist(lapply(GSE55389.magenta.filter,ROC_plot2,data=GSE55389.data,y="typeII"))
  dev.off()
  GSE55389.magenta.ROC=data.frame(gene=GSE55389.magenta.filter,GSE55389.ROC=GSE55389.magenta.ROC)
  magenta.ROC=merge(magenta.filter,GSE55389.magenta.ROC,by = "gene",all.x = T)
  write.csv(magenta.ROC,"magenta.ROC.csv")

}


# boxplot for Itgb3bp
datExpr[1:4,1:4]
exprSet = t(datExpr)
head(datTraits)
design=model.matrix(~0+ datTraits$type)
design=data.frame(design,design[,2]+design[,3])
colnames(design)=c(levels(factor(datTraits$type)),"both_type")

usedata=data.frame(datExpr,design)

gene_box<-function(gene="Tmbim1",group="both_type",data=usedata){
  p <- ggboxplot(data, x = group, y = gene,
                 ylab = sprintf("Expression of %s",gene),
                 xlab = group,
                 color = group, 
                 palette = "nejm",
                 add = "jitter")
  p + stat_compare_means(label = "p.format",label.x=1.5)
}
A=gene_box(gene="Tmbim1",group="both_type",data=usedata)
B=gene_box(gene="Col4a2",group="typeII",data=usedata)
C=gene_box(gene="Cdan1",group="typeII",data=usedata)

pdf("Itgb3bp_in_traindataset.pdf", useDingbats=FALSE,width = 2.5,height = 2.5)
gene_box(gene="Itgb3bp ",group="typeII",data=usedata)
dev.off()
