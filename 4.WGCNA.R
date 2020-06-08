rm(list = ls())
options(stringsAsFactors = F)



library('gplots')
library('ggplot2')
library('knitr')
library('reshape2')
library('RColorBrewer')
library('WGCNA')

###WGCNA

{
  load("filtered.Rdata")
  enableWGCNAThreads()
  rownames(exprs.upper)=gene_map$SYMBOL[match(rownames(exprs.upper),gene_map$ENTREZID)]
  datExpr=t(exprs.upper)
  datExpr[1:4,1:4]
  datTraits = data.frame(sample = rownames(pdata),pdata)
  
  ## step 1 :
  if(T){
    datExpr[1:4,1:4]
    head(datTraits)
    dim(datExpr)
  
    gsg = goodSamplesGenes(datExpr, verbose = 3)
    
    if (!gsg$allOK){
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes)>0) 
        printFlush(paste("Removing genes:", 
                         paste(names(WGCNA_matrix)[!gsg$goodGenes], collapse = ",")));
      if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", 
                         paste(rownames(WGCNA_matrix)[!gsg$goodSamples], collapse = ",")));
      # Remove the offending genes and samples from the data:
      WGCNA_matrix = WGCNA_matrix[gsg$goodSamples, gsg$goodGenes]
    }
  }
  
  #step 2
  # remove outliers
  if(T){
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    #首先针对样本做个系统聚类
    datExpr_tree<-hclust(dist(datExpr), method = "average")
    par(mar = c(0,5,2,0))
    plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
         cex.axis = 1, cex.main = 1,cex.lab=1)
    # threshold 25
    abline(h = 25, col = "red") 
    clust = cutreeStatic(datExpr_tree, cutHeight = 25, minSize = 10)
    table(clust) # 0 means cut，1 means keep
    keepSamples = (clust==1|2)
    datTraits = datTraits[keepSamples, ]
    datExpr = datExpr[keepSamples,]
    sampleTree2 = hclust(dist(datExpr), method = "average")
    plot(sampleTree2, main = "Sample clustering to detect outliers", 
         sub="", xlab="", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2)
    
    ## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
    #针对前面构造的样品矩阵添加对应颜色
    type_colors <- numbers2colors(as.numeric(factor(datTraits$type)), 
                                  colors = c("grey","green","red"),signed = F)
    disease_colors <- numbers2colors(as.numeric(factor(datTraits$disease)), 
                                     colors = c("green","grey","red","white"),signed = F)
    ## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目
    # age_colors <- numbers2colors( as.numeric(datTraits$age) ,signed = FALSE)
    ## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，当然，这样给的颜色不明显，意义不大
    #10个样品的系统聚类树及性状热图
    merge_color = cbind(type_colors)
    dimnames(merge_color) = dimnames(pdata[,6])
    par(mar = c(1,4,3,1),cex=0.8)
    
    
    pdf("sample-subtype-cluster.pdf",width = 8,height = 6,useDingbats = F)
    plotDendroAndColors(datExpr_tree, merge_color,
                        groupLabels = colnames(sample),
                        cex.dendroLabels = 0.8,
                        marAll = c(1, 4, 3, 1),
                        cex.rowText = 0.01,
                        main = "Sample dendrogram and trait heatmap")
    # abline(h=25, col = "red")
    dev.off()
  }
  
  
  
  ## step 3   ##determine power
  if(T){
    powers = c(seq(1, 10, by=1),seq(12, 20, by=2))
    # Call the network topology analysis function
    sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    sft$powerEstimate
    # [1] 7
    
    png("step2-beta-value.png",width = 1200,height = 600)
    # Plot the results:
    ##sizeGrWindow(9, 5)
    par(mfrow = c(1,2));
    cex1 = 0.9;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    ##横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
    ##网络越符合无标度特征 (non-scale)
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.9,col="red")
    ##筛选标准。R-square=0.9
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
    
    # power=7
    
    
    ADJ1_cor <- abs(WGCNA::cor(datExpr,use = "p" ))^7
    k <- as.vector(apply(ADJ1_cor,2,sum,na.rm=T))
    # 基因多的时候使用下面的代码：
    # k <- softConnectivity(datE=datExpr,power=7) 
    sizeGrWindow(10, 5)
    pdf("power_second_graph.pdf",width = 12,height = 6, useDingbats=FALSE)
    par(mfrow=c(1,2))
    hist(k,xlab = "Connectivity")
    scaleFreePlot(k,main="Check Scale free topology\n")
    dev.off()
  }
  
  ## step4 construction of co-expression network
  if(T){
    #(1)Co-expression similarity and adjacency 
    adjacency = adjacency(datExpr, power = 7) 
    #(2) Turn adjacency into topological overlap
    TOM = TOMsimilarity(adjacency);
    dissTOM = 1-TOM
    # (3)  Call the hierarchical clustering function
    geneTree = hclust(as.dist(dissTOM), method = "average");
    # Plot the resulting clustering tree (dendrogram)
    sizeGrWindow(12,9)

    plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
         labels = FALSE, hang = 0.04);
    # We like large modules, so we set the minimum module size relatively high:
    minModuleSize = 30;
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize);
    table(dynamicMods)
    # Convert numeric lables into colors
    dynamicColors = labels2colors(dynamicMods)
    table(dynamicColors)
    # Plot the dendrogram and colors underneath
    #sizeGrWindow(8,6)
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    
    # Calculate eigengenes
    MEList = moduleEigengenes(datExpr, colors = dynamicColors)
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs);
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average");
    # Plot the result
    #sizeGrWindow(7, 6)
    plot(METree, main = "Clustering of module eigengenes",
         xlab = "", sub = "")
    #merge modules with 75% similarity
    MEDissThres = 0.25
    # Plot the cut line into the dendrogram
    abline(h=MEDissThres, col = "red")
    # Call an automatic merging function
    merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors = merge$colors;
    table(mergedColors)
    # Eigengenes of the new merged modules:
    mergedMEs = merge$newMEs
    pdf("Gene dendrogram and module colors.pdf",width = 8,height = 6, useDingbats=FALSE)
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    dev.off()
    net=merge
    save(net,datExpr,datTraits,gene_map,file="WGCNA.Rdata")
  }

  load("WGCNA.Rdata")
  
  
  
  
  ## step 5 
  ## modul-traits heatmap
  if(T){
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    design=model.matrix(~0+ datTraits$type)
    colnames(design)=c(levels(factor(datTraits$type)))
     moduleColors <- net$colors##extract colors
    # Recalculate MEs with color labels
    moduleEG= moduleEigengenes(datExpr, moduleColors)
    MEs0 = moduleEG$eigengenes
    MEs = orderMEs(MEs0)
    plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), plotDendrograms = T, xLabelsAngle = 90)
    # 根据基因间表达量进行聚类所得到的各模块间的相关性图
    # marDendro/marHeatmap 设置下、左、上、右的边距
    moduleTraitCor = cor(MEs, design , use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    
    sizeGrWindow(10,10)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    pdf("step5-12Module-trait-relationships.pdf",width = 8,height = 12, useDingbats=FALSE)
    par(mar = c(6, 8.5, 3, 3));
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(design),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.8,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    dev.off()
    
    # identification of most significant moduls
    if (T) {
      cor_ADR <- signif(WGCNA::cor(design,MEs,use="p",method="pearson"),5)
      names(cor_ADR)= rep(colnames(cor_ADR),each=nrow(cor_ADR))
      p.values <- corPvalueStudent(cor_ADR,nSamples=nSamples)
      names(p.values)= rep(colnames(p.values),each=nrow(p.values))
      Freq_MS_max_cor <- names(cor_ADR[which.max(abs(cor_ADR[,-which(colnames(cor_ADR) == "MEgrey")]))])
      Freq_MS_max_cor
      Freq_MS_max_p <- names(p.values[which.min(p.values[,-which(colnames(p.values) == "MEgrey")])])
      Freq_MS_max_p
    }
    
  
  
  # step6 
  # TOM plot, take a long time
  if(F){
        ## Loading objects:
    ##   TOM
    
    TOM <- as.matrix(TOM)
    
    dissTOM = 1-TOM
    
    selectTOM = dissTOM
    selectTree = hclust(as.dist(selectTOM), method = "average")
    selectColors = moduleColors;
    # Open a graphical window
    sizeGrWindow(9,9)
    # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
    # the color palette; setting the diagonal to NA also improves the clarity of the plot
    plotDiss = selectTOM^7;
    diag(plotDiss) = NA;
    myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
    png("step7-All-Network-heatmap.png",width = 800,height = 600)
    TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",col=myheatcol)
    dev.off()
  }
  
  
  # step 7
  # gene counts plot & extrat genes in interested modules
  if(T){
    ##
    pdf("gene-counts.pdf",width = 8,height = 6, useDingbats=FALSE)
    x = barplot(table(moduleColors),ylab = "gene counts",xlab = "module",
                col =names(table(moduleColors)),ylim = c(0,2000),
                main = "gene counts in all modules")
    text(x,table(moduleColors),labels  = table(moduleColors),cex=1.5,pos=3)
    dev.off()
    seleteced.module.gene=function(x){
      module = x
      # Select module probes
      probes = colnames(datExpr) ## 
      inModule = (moduleColors==module);
      modProbes = probes[inModule];  
      return(modProbes)
    }
    related.module = c("turquoise","blue","magenta","greenyellow",
                       "pink","purple")
    for (i in related.module) {
      module = i
      gene = seleteced.module.gene(module)
      write.table(gene,paste(module,"selected-module-gene.txt",sep = "-"),row.names=F,col.names = F,quote=F)
    }
    for (i in related.module) {
      module = i
      gene = seleteced.module.gene(module)
      gene.entriz=gene_map$ENTREZID[match(gene,gene_map$SYMBOL)]
      # tangene  = seleteced.module.gene("tan")
      # greengene = seleteced.module.gene("green")
      # save(bluegene,file = "selected-module-gene.Rdata")
      write.table(gene.entriz,paste(module,"selected-module-ENTREZID.txt",sep = "-"),row.names=F,col.names = F,quote=F)
    }
  }
    
  ## step 8 
  # export connection to cytoscape
  if(T){
    # Recalculate topological overlap
    TOM = TOMsimilarityFromExpr(datExpr, power = 7); 
    related.module = c("blue","magenta","turquoise")
    for (i in related.module) {
      # Select module
      module = i
      # Select module probes
      probes = colnames(datExpr) 
      inModule = (moduleColors==module);
      modProbes = probes[inModule]; 
      # Select the corresponding Topological Overlap
      modTOM = TOM[inModule, inModule];
      dimnames(modTOM) = list(modProbes, modProbes)
      ## 
      cyt = exportNetworkToCytoscape(
        modTOM,
        edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
        nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
        weighted = TRUE,
        threshold = 0.02,
        nodeNames = modProbes, 
        nodeAttr = moduleColors[inModule]
      );
    }
    
  }
}



