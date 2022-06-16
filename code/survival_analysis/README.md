#### Survival analysis

Author: [Komal Rathi](https://github.com/komalsrathi)

Code: [survival_analysis](../survival_analysis)

##### Introduction 

The goal of this module is to compare the survival of samples that cluster near vs those that cluster away from the patient of interest using the R package `survival`. 

#####  Scripts

1. `p5_kaplan_meier_pediatric.R` and `p6_kaplan_meier_adult.R`: Using the distance metric from umap, two clusters are determined: the top 20 nearest samples from umap analysis are treated as one cluster i.e. `Cluster With Patient` and all other comparator tumors are treated as a second cluster i.e. `Cluster away from Patient`. Using this information, a kaplan-meier plot showing the overall survival between the two groups is plotted.