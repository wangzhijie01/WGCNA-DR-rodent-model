rm(list = ls())
options(stringsAsFactors = F)

library(GEOquery)
library(stringr)
getGEOdata <- function(x, raw=FALSE){
  # Use the getGEO() function to download the GEO data for the id stored in x
  GSEdata <- getGEO(x, GSEMatrix=T, AnnotGPL=FALSE)
  # check the data: printing a summary of the expression values for the first 2 columns
  print(summary(exprs(GSEdata[[1]])[,1:2]))
  # Get the eset object
  eset <- GSEdata[[1]]
  # Save the objects generated for future use in the current working directory
  save(GSEdata, eset, file=paste(x, ".RData", sep=""))
  # check whether we want to return the list object we downloaded on GEO or
  # just the eset object with the getGSEobject argument
  if(raw) return(GSEdata) else return(eset)
}
GSE111465=getGEOdata("GSE111465")
GSE28831=getGEOdata("GSE28831")
GSE55389=getGEOdata("GSE55389")
save(GSE111465,GSE28831,GSE55389,file="raw.GSE.Rdata")

load("raw.GSE.Rdata")
myGSE=data.frame(GSE.ID=c("GSE111465","GSE28831","GSE55389"),
                 GPL.ID=c("GPL20710","GPL7294","GPL6246"))

#convet to SYMBOL
{
# ID convert function
useGPLtoID=function(GSE.ID=GSE.ID,GPL.ID=GPL.ID,IDname=IDname,SYMBOLname=SYMBOLname)
  {
    load(file = paste(GSE.ID,".Rdata",sep=""))
    exprSet=exprs(eset)
    range(exprSet)
    exprSet=as.data.frame(exprSet)
    #用GPL文件注释
    library(GEOquery)
    gpl = getGEO(GPL.ID, destdir="./Annotation files/")
    names(Meta(gpl))
    colnames(Table(gpl))
    colname=str_subset(colnames(Table(gpl)),pattern = ".*(Symbol|symbol|SYMBOL).*")[1]
    iddata <- Table(gpl)[c("ID",colname)]
    head(iddata)
    mann.anno=iddata[match(rownames(exprSet),iddata[,"ID"]),]
    exprSet$Symbol=mann.anno[,colname]
    exprSet=na.omit(exprSet)
    dat <- aggregate(x = exprSet[,1:(ncol(exprSet)-1)],
                     by = list(exprSet$Symbol),
                     FUN = mean)
    rownames(dat) = dat$Group.1
    dat=dat[,-1]
    dim(dat)
    dat[1:4,1:4]
    save(dat,file=paste(GSE.ID,"filtered.Rdata",sep = "."))
    return(dat)
    print(paste(GSE.ID,"finished",sep = " "))
  }

  # "GSE111465","GSE28831",
allGSE <- lapply(1:2,
                 function(i)
                 {
                   GSE.ID=myGSE[i,1]
                   GPL.ID=myGSE[i,2]
                   temp_list = useGPLtoID(GSE.ID,GPL.ID)
                   return(temp_list)
                 })
# "GSE55389"
{
    load(file = paste("GSE55389",".Rdata",sep=""))
    exprSet=exprs(eset)
    range(exprSet)
    exprSet=as.data.frame(exprSet)
   
    library(GEOquery)
    gpl = getGEO("GPL6246", destdir="./Annotation files/")
    names(Meta(gpl))
    colnames(Table(gpl))
    gpl_anno  <- Table(gpl)
    library(dplyr)
    library(tidyr)
    probe2symbol_df <- gpl_anno %>%
      
      dplyr::select(ID,gene_assignment) %>%
      
      filter(gene_assignment != '---') %>%
      
      separate(gene_assignment,c('drop','symbol'),sep=' // ') %>%
      
      dplyr::select(-drop)
    
    head(probe2symbol_df)
    
    mann.anno=probe2symbol_df[match(as.numeric(rownames(exprSet)),probe2symbol_df[,"ID"]),]
    exprSet$Symbol=mann.anno[,"symbol"]
    exprSet=na.omit(exprSet)
    dat <- aggregate(x = exprSet[,1:(ncol(exprSet)-1)],
                     by = list(exprSet$Symbol),
                     FUN = mean)
    rownames(dat) = dat$Group.1
    dat=dat[,-1]
    dim(dat)
    dat[1:4,1:4]
    save(dat,file=paste("GSE55389","filtered.Rdata",sep = "."))
    print(paste("GSE55389","finished",sep = " "))
}
allGSE[[3]]=dat
# extract common symbols
load("filtered.Rdata")
colno=match(rownames(exprs.upper),gene_map$ENTREZID)
ref.symbol=gene_map$SYMBOL[colno]
for (i in 1:length(allGSE)) {
  tmp=allGSE[[i]]
  col=match(ref.symbol,rownames(tmp),nomatch =0 )
  allGSE[[i]]=allGSE[[i]][col,]
}
}
names(allGSE)=myGSE$GSE.ID
save(allGSE,file = "allGSE.Rdata")


