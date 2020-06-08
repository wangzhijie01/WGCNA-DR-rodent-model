seleteced.module.gene=function(x){
  module = x
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  inModule = (moduleColors==module);
  modProbes = probes[inModule];  
  return(modProbes)
}