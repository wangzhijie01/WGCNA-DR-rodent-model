rm(list = ls())
options(stringsAsFactors = F)


load(file = "aftercombat.Rdata")
exprSet=as.data.frame(exprSet2)
range(exprSet)
dev.off()
boxplot(exprSet,outline = F) #check 

#ID convert
if (T) {
  library(GEOquery)
  gpl = getGEO('GPL6247', destdir="./Annotation files/")
  names(Meta(gpl))
  colnames(Table(gpl))
  # 选取其中的几列用来转换
  gpl_anno  <- Table(gpl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  iddata <- gpl_anno %>%
    dplyr::select(ID,gene_assignment) %>%
    filter(gene_assignment != '---')

  tmp=str_match(iddata$gene_assignment,pattern = '(\\w+) // (\\w+) // (.*) // ([0-9]+$)')
  
  iddata=data.frame(ID=iddata$ID,Symbol=tmp[,3],ENTREZID=tmp[,5])
  
  iddata$ID=as.character(iddata$ID)
  head(iddata)
  mann.anno=iddata[match(rownames(exprSet),iddata$ID),]
  exprSet$ENTREZID=mann.anno$ENTREZID
  exprSet=na.omit(exprSet)
  # mean value
  dat <- aggregate(x = exprSet[,1:(ncol(exprSet)-1)],
                   by = list(exprSet$ENTREZID),
                   FUN = mean)
  rownames(dat) = dat$Group.1
  dat=dat[,-1]
  dim(dat)
  exprSet=dat
}

#filter no ENTREZID gene
library('org.Rn.eg.db')
gene_map<-mapIds(org.Rn.eg.db, as.character(rownames(exprSet)), 'SYMBOL', keytype ='ENTREZID',multiVals = "first")
gene_map=data.frame(ENTREZID=names(gene_map),SYMBOL=as.character(gene_map))
gene_map=na.omit(gene_map)
exprSet=exprSet[rownames(exprSet)%in%gene_map$ENTREZID,]

table(is.na(gene_map$SYMBOL))
{
  #get top 25% MAD genes
  m.mad=apply(exprSet,1,mad)
  exprs.upper=exprSet[which(m.mad>quantile(m.mad, probs = seq(0, 1, 0.25))[4]),]
  exprs.upper=as.data.frame(exprs.upper)
  
}
save(exprs.upper,exprSet,pdata,gene_map,file="filtered.Rdata")

