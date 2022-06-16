#### Immune analysis

Author: [Komal Rathi](https://github.com/komalsrathi)

Code: [immune_analysis](../immune_analysis)

##### Introduction 

The goal of this module is to perform immune deconvolution on input RNA-seq expression matrices using the R package `xCell`. 

#####  Scripts

1. `p4_immune_profile_pediatric.R` and `p4_immune_profile_adult.R`: All protein codings genes from the batch corrected and back-transformed count matrix are used for immune deconvolution using `xCell`. Next, a boxplot representing median immune-scores per cell type is plotted for all comparator tumors (either pediatric or adult tumors) and the patient of interest, with the latter highlighted in red.

2. `p5_immune_profile_topcor_pediatric.R`: All protein codings genes from the batch corrected and back-transformed count matrix are used for immune deconvolution using `xCell`. Next, a boxplot representing median immune-scores per cell type is plotted for the top 20 similar pediatric tumors obtained from umap and the patient of interest, with the latter highlighted in red.

3. `p4_tis_profile.R`: Here, a combined counts matrix filtered to only protein coding genes from all pediatric tumors, adult tumors and patient of interest is created. To reduce the size of the matrix only cancer groups with >=5 samples are kept. A quantile normalization is performed after log-transforming the merged count matrix. Next, the normalized matrix is filtered to an 18-gene tumor inflammation signature (TIS) score to provide a readout of key immune functions (cytoxicity, exhaustion, immune checkpoint expression, etc). This signature is obtained from: https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-019-2100-3#Sec14. Finally, a boxplot representing average tumor inflammation score per cancer group is plotted with the y-intercept corresponding to the tumor inflammation score of the patient of interest demarcating where the POI falls relative to adult and pediatric tumors.

