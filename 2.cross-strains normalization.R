rm(list = ls())
options(stringsAsFactors = F)


library(sva)
library('gplots')
library('ggplot2')
library('knitr')
library('reshape2')
library('RColorBrewer')
library('WGCNA')


#####extract T1DM OR T2DM
if(T){
  load(file = 'eset.Rdata')
  
  rownames(pdata)
  pdata = pdata[1:23,]
  exprSet = exprSet[,1:23]
}

####check,PCA
if(T){
  group_list = pdata$strain
  dat = exprSet
  table(group_list)
  dat[1:4,1:4]
  dat=t(dat)
  dat=as.data.frame(dat)
  dat=cbind(dat,group_list) 
  library("FactoMineR")
  library("factoextra") 
  dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = dat$group_list, # color by groups
               #palette = c("#00AFBB", "#E7B800"),
               addEllipses = F, # Concentration ellipses
               legend.title = "strain"
  )
  ggsave('all_samples_PCA.pdf',width = 4,height = 4, useDingbats=FALSE)
}


##cross-strains normalization
{
  modcombat = model.matrix(~as.factor(type), data=pdata)#
  batch1 = pdata[,1]
  combat_edata = ComBat(dat = as.matrix(exprSet), batch=batch1, par.prior=TRUE, prior.plots=F)
  exprSet2 = combat_edata
  # hcluster
  sampleTree2 = hclust(dist(t(exprSet2)), method = "average")
  png("all_samples_cluster_after_combact.png",width = 800,height = 600)
  plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="")
  dev.off()
  
  group_list = pdata$strain
  dat = exprSet2
  table(group_list)
  dat[1:4,1:4]
  dat=t(dat)#
  dat=as.data.frame(dat)
  dat=cbind(dat,group_list)
  library("FactoMineR")
  library("factoextra") 
  dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = dat$group_list, # color by groups
               #palette = c("#00AFBB", "#E7B800"),
               addEllipses = F, # Concentration ellipses
               legend.title = "strain")
  ggsave('all_samples_PCA_after_combact.pdf',width = 4,height = 4, useDingbats=FALSE)
  save(pdata,exprSet2,file = "aftercombat.Rdata")
}


