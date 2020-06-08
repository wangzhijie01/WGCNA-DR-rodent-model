rm(list = ls())
options(stringsAsFactors = F)
# Load the package
library(WGCNA)
load("WGCNA.Rdata")
source("select.module.gene.R")

# intramodularConnectivity chose hub gene
  if(T){
    adjacency = adjacency(datExpr, power = 7) 
    moduleColors=net$colors
    # Convert labels to colors for plotting
    # mergedColors = labels2colors(net$colors)
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    design=model.matrix(~0+ datTraits$type)
    colnames(design)=levels(factor(datTraits$type))
    # Recalculate MEs with color labels
    moduleEG= moduleEigengenes(datExpr, moduleColors)
    MEs0 = moduleEG$eigengenes
    
    MEs = orderMEs(MEs0); 
    MM= as.data.frame(cor(datExpr, MEs, use ="p"))
    KIM = intramodularConnectivity(adjacency, moduleColors, scaleByMax= TRUE)
    
    datKME=signedKME(datExpr, MEs, outputColumnName="kME_MM.")	
    
    HubGenes <- chooseTopHubInEachModule(datExpr,moduleColors)
    write.table(HubGenes,file = "HubGenes_of_each_module.txt",quote=F,sep='\t')
  }

# hub genes plots
{
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  design=model.matrix(~0+ datTraits$type)
  colnames(design)=c(levels(factor(datTraits$type)))
  moduleColors=net$colors
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  hubgene=c("Col4a2","Cdan1","Tmbim1")
  x=geneModuleMembership[hubgene,]
  type = as.data.frame(design[,2:3]);
  type$both=type[,1]+type[,2]
  
  
  geneTraitSignificance = as.data.frame(cor(datExpr, type, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  # names(geneTraitSignificance) = paste("GS.", names(type), sep="");
  # names(GSPvalue) = paste("p.GS.", names(type), sep="");
  
  related.module = c("blue","magenta","turquoise")
  related.trait=c("both","typeII","typeII")
  scater.data=data.frame()
  for (i in 1:3) {
    module = related.module[i]
    column = match(module, modNames);
    moduleGenes =hubgene[i]
    trait=related.trait[i]
    scater.data[i,"MM"]=abs(geneModuleMembership[moduleGenes, column])
    scater.data[i,"GS"]=abs(geneTraitSignificance[moduleGenes,trait])
    scater.data[i,"GS.p"]=GSPvalue[moduleGenes,trait]
  }
  rownames(scater.data)=hubgene
  scater.data$col=related.module
  scater.data$ID= rownames(scater.data)
  pdf(paste("typeII-related",module,"step6-Module_membership-gene_significance.pdf",sep = "-"),width = 6,height = 4,useDingbats = F)
  p=ggplot(data = scater.data,aes(x = MM, y = GS, colour = col)) + 
    geom_point(size = 3.5)+
    scale_color_manual(values = scater.data$col)+
    xlim(0.5,1)+
    ylim(0,0.7)+
    labs(x = "Module Membership", y = "Gene significance")+
    ggtitle("Module membership vs. gene significance")+theme_bw()
  
  p +
    geom_point(size = 3, shape = 1, data = scater.data) +
    ggrepel::geom_label_repel(
      aes(label = ID),
      data = scater.data,
      color="black"
    )
  
  dev.off()
}
