rm(list = ls())
options(stringsAsFactors = F)
# Load the package
library(WGCNA);

#Read in the female liver data set
# The following 3421 probe set were arrived at using the following steps
#1) reduce to the 8000 most varying, 2) 3600 most connected, 3) focus on unique genes
#reference dataset
load("WGCNA.Rdata")
datExpr[1:4,1:4]
head(gene_map)
colno=match(colnames(datExpr),gene_map$ENTREZID)
colnames(datExpr)=gene_map$SYMBOL[colno]
datExpr.ref=datExpr

# color.ref = labels2colors(net$colors)
color.ref = net$colors


#compare datasets
load("allGSE.Rdata")
mp.all=lapply(1:length(allGSE), function(i) {
  datExpr.compare=allGSE[[i]]
  datExpr.compare=as.data.frame(t(datExpr.compare))
  datExpr.compare[1:4,1:4]
  setLabels = c("ref", "compare");
  multiExpr = list(ref = list(data = datExpr.ref), compare = list(data = datExpr.compare));
  multiColor = list(ref = color.ref);
  mp = modulePreservation(multiExpr, multiColor,
                                      referenceNetworks = 1,
                                      nPermutations = 200,
                                      randomSeed = 1,
                                      quickCor = 0,
                                      verbose = 3)
  
  return(mp)
})

  names(mp.all)=names(allGSE)
  save(mp.all,file="mp.all.Rdata")


  load(file = "mp.all.Rdata")
  for (i in 1:length(mp.all)) {
    mp=mp.all[[i]]
    setLabels = c("ref", "compare");
    ref = 1
    test = 2
    statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
    statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
    # Compare preservation to quality:
    print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                 signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
    
    # Module labels and module sizes are also contained in the results
    modColors = rownames(mp$preservation$observed[[ref]][[test]])
    moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
    # leave grey and gold modules out
    plotMods = !(modColors %in% c("grey", "gold"));
    # Text labels for points
    text = modColors[plotMods];
    # Auxiliary convenience variable
    plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
    # Main titles for the plot
    mains = c("Preservation Median rank", "Preservation Zsummary");
    # Start the plot
    sizeGrWindow(10, 5);
    filename=paste(names(mp.all)[i],"BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf")
    pdf(file=filename, wi=10, h=5, useDingbats=FALSE)
    par(mfrow = c(1,2))
    par(mar = c(4.5,4.5,2.5,1))
    for (p in 1:2)
    {
      min = min(plotData[, p], na.rm = TRUE);
      max = max(plotData[, p], na.rm = TRUE);
      # Adjust ploting ranges appropriately
      if (p==2)
      {
        if (min > -max/10) min = -max/10
        ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
      } else
        ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
      plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
           main = mains[p],
           cex = 2.4,
           ylab = mains[p], xlab = "Module size", log = "x",
           ylim = ylim,
           xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
      labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08,protectEdges=F);
      # For Zsummary, add threshold lines
      if (p==2)
      {
        abline(h=0)
        abline(h=2, col = "blue", lty = 2)
        abline(h=10, col = "darkgreen", lty = 2)
      }
      }
      # If plotting into a file, close it
      dev.off();
    
    } 

# Save the results





# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
sizeGrWindow(12, 9);
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
# Plot each Z statistic in a separate plot.
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/5) min = -max/5
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.7,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(20, 1000))
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}

data.frame(color = modColors[plotMods], label = labs)