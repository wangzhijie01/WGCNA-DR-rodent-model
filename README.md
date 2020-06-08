# WGCNA-DR-rodent-model
------
> This is the raw code for paper "The characteristic Networks based on Weighted Gene Co-expression Network Analysis in Diabetic Retinopathy Rodent Models".

Raw data was download from [E-MTAB-5563](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5563/)

Other three datasets for Module Preservation Analysis were download from [GSE28831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28831), [GSE111465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111465), [GSE55389](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55389).

For repeating results, you should run the code for the following order.
- 1.affyexpress.R
- 2.cross-strains normalization.R
- 3.gene filter.R
- 4.WGCNA.R
- 5.other datasets.R
- 6.modulePreservation.R
- 7.hub-genes.R
- 8.ROC.R
