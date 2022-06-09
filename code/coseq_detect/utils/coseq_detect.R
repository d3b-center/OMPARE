# Author: Run Jin
# Run COSEQ

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(coseq)
  library(diptest)
  library(pheatmap)
  library(NB.MClust)
  library(edgeR)
  library(rlist)
  library(RColorBrewer)
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("--ref_cancer_dir"),type="character",
              help="directory with cancer type specific reference tumor subset"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory of patient of interest")
)
opt <- parse_args(OptionParser(option_list=option_list))
ref_cancer_dir <- opt$ref_cancer_dir
output_dir <- opt$output_dir

# Read in files necessary for analyses
histology_df <- list.files(path = ref_cancer_dir, pattern = "histologies", full.names = T)
histology_df <- readr::read_tsv(histology_df)
exp_count_cg_coding <- list.files(path = ref_cancer_dir, pattern = "count", full.names = T)
exp_count_cg_coding <- readRDS(exp_count_cg_coding)

# give them the same group for all samples
group <- factor(rep(c(1), times = length(colnames(exp_count_cg_coding))))

# filter lowly expressed genes by DESeq2
exp_count_cg_coding <- DGEList(counts = exp_count_cg_coding, group = group)
keep <- filterByExpr(exp_count_cg_coding)
exp_count_cg_coding <- exp_count_cg_coding[keep,,keep.lib.sizes=FALSE] 

# extract the resulting counts from the DGElist object as a df to use for dip.test
filtered_counts_cg_coding <- as.data.frame(exp_count_cg_coding$counts)
saveRDS(filtered_counts_cg_coding, file = file.path(output_dir, "expression_matrix.rds"))

# compute the normalization factors (TMM) for scaling by library size
exp_count_cg_coding.tmm <- calcNormFactors(exp_count_cg_coding)

#### DIP TEST
# define a dataframe to store the pval
filtered_count_cg_coding_pval <- data.frame()
for(i in 1:nrow(filtered_counts_cg_coding)){
  # calculate dip.test for each gene
  diptest_each <- filtered_counts_cg_coding[i,] %>% 
    as.numeric() %>%
    dip.test(simulate.p.value = FALSE, B = 2000)
  # gather the pvalue for each gene
  pval_each <- diptest_each$p.value
  # add another column to store pval
  filtered_count_cg_pval_each <- filtered_counts_cg_coding[i,] %>%
    mutate(pval = pval_each)
  # combine each line back to the dataframe 
  filtered_count_cg_coding_pval <- bind_rows(filtered_count_cg_pval_each,  filtered_count_cg_coding_pval)
}

# Filter to expression that only has pval<0.05 if over 1000 genes satisfy this or take top 1000
# see how many genes have pval<0.5
n <- filtered_count_cg_coding_pval %>% 
  dplyr::filter(pval <0.05) %>%
  nrow()

if(n>1000){
  filtered_count_cg_coding_pval <-filtered_count_cg_coding_pval %>% 
    dplyr::filter(pval <0.05) %>% 
    dplyr::select(-pval)
} else {
  filtered_count_cg_coding_pval <-filtered_count_cg_coding_pval %>% 
    dplyr::arrange(pval, descending = FALSE) %>% 
    head(1000) %>%
    dplyr::select(-pval)
}

# generate DGEList using the pval, protein coding and gene expression filtered count
filtered_cpm_cg_coding_pval <- DGEList(counts = filtered_count_cg_coding_pval, 
                                       # norm factor calculated on the entire dataset
                                       norm.factors = exp_count_cg_coding.tmm$samples$norm.factors, 
                                       # library size calculated on the entire dataset
                                       lib.size = colSums(exp_count_cg_coding.tmm$counts))
# calculate CPM and use the CPM as input for COSEQ analysis 
filtered_cpm_cg_coding_pval <-cpm(filtered_cpm_cg_coding_pval, log = FALSE) %>%
  as.data.frame()

# Run COSEQ on dip.test + CPM results with Logit transformation
runLogit <- coseq(filtered_cpm_cg_coding_pval, 
                  K=2:20, 
                  model="Normal", 
                  transformation="logit", 
                  GaussianModel = "Gaussian_pk_Lk_I",
                  seed=2021)

# save the plots for diagnosis
pdf(file.path(output_dir, "logit_Gaussian_pk_Lk_I_plot.pdf"))
coseq::plot(runLogit)
dev.off()

# save the results 
sink(file.path(output_dir, "logit_Gaussian_pk_Lk_I_summary.txt"), type=c("output"))
summary(runLogit) 
sink()

# obtain cluster information
cluster_coseq <- coseq::clusters(runLogit) %>% 
  as.data.frame() 

############################### utilize NB.Mclust to generate sample cluster
# convert to numeric and round
filtered_count_cg_coding_pval[] <- lapply(filtered_count_cg_coding_pval, function(x) {
  if(is.factor(x)){
    round(as.numeric(as.character(x)))
  }  else {
    round(as.numeric(x))
  }
})

# save output
saveRDS(filtered_count_cg_coding_pval, file = file.path(output_dir, "filtered_count_cg_coding_pval.rds"))

# transpose the matrix for fitting the data 
filtered_tcount_cg_coding_pval <- t(filtered_count_cg_coding_pval)

# run NB.mClust
nb_mclust_res <- NB.MClust(filtered_tcount_cg_coding_pval,
                           K=2:10, 
                           tau0=10,
                           rate=0.9,
                           bic=TRUE,
                           iteration=50)
# save optimal K 
sink(file.path(output_dir, "nb_clust_k_selected.txt"), type=c("output"))
nb_mclust_res$K 
sink()

############################### plot heatmap with NB.MClust clustering results
#Estimated cluster assignment 
cluster_nb <- nb_mclust_res$cluster %>%
  as.data.frame()
colnames(cluster_nb) <- c("cluster_assigned_nb")

####### take log and then zscore of the matrix for better visualization
filtered_log_cg_coding_pval<-log2(filtered_count_cg_coding_pval+1)
# calculate mean
matrix_means <-rowMeans(filtered_log_cg_coding_pval, na.rm = TRUE)
# calculate sd
matrix_sd <- apply(filtered_log_cg_coding_pval, 1, sd, na.rm = TRUE)
# subtract mean
filtered_log_cg_coding_pval_minus_mean <- sweep(filtered_log_cg_coding_pval, 1, matrix_means, FUN = "-")
# divide by SD remove NAs and Inf values from zscore for genes with 0 in normData
filtered_log_cg_coding_pval_zscored <- sweep(filtered_log_cg_coding_pval_minus_mean, 1,matrix_sd, FUN = "/") %>% 
  na_if(Inf) %>% na.omit()

############## arrange annotation lists by cluster and molecular subtype information
# arrange coseq annotation file
colnames(cluster_coseq) <- c("cluster_by_coseq")
cluster_coseq <- cluster_coseq %>%
  dplyr::arrange(cluster_by_coseq) %>%
  tibble::rownames_to_column("gene_symbol")
cluster_coseq$cluster_by_coseq <- as.factor(cluster_coseq$cluster_by_coseq) 

# arrange the genes based on gene cluster by COSEQ
filtered_log_cg_coding_pval_zscored <- filtered_log_cg_coding_pval_zscored %>%
  tibble::rownames_to_column("gene_symbol") %>% 
  left_join(cluster_coseq) %>% 
  arrange(cluster_by_coseq) %>% 
  dplyr::select(-cluster_by_coseq) %>% 
  tibble::column_to_rownames("gene_symbol")

# arrange NB.MClust cluster 
cluster_nb <- cluster_nb %>%
  dplyr::arrange(cluster_assigned_nb) %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID")
cluster_nb$cluster_assigned_nb <- as.factor(cluster_nb$cluster_assigned_nb)

# generate subtype annotation file annotation file
subtype_anno <- histology_df %>% dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype) %>% 
  mutate(molecular_subtype = as.character(molecular_subtype)) %>%
  mutate(molecular_subtype = case_when(
    is.na(molecular_subtype) ~ "Not available",
    TRUE ~ molecular_subtype
  ))

# combine subtype and cluster from NB.Mclust order
combined_anno <- cluster_nb %>% 
  left_join(subtype_anno) %>%
  arrange(molecular_subtype) %>%
  arrange(cluster_assigned_nb)
combined_anno$molecular_subtype <- as.factor(combined_anno$molecular_subtype)
combined_anno$cluster_assigned_nb <- as.factor(combined_anno$cluster_assigned_nb)

# arrange the matrix based on molecular subtype cluster assigned by NB.MClust
filtered_log_cg_coding_pval_zscored_anno <- filtered_log_cg_coding_pval_zscored %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::left_join(combined_anno) %>%
  arrange(molecular_subtype) %>%
  arrange(cluster_assigned_nb) %>%
  dplyr::select(-molecular_subtype) %>%
  dplyr::select(-cluster_assigned_nb) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")


# revert back annotation columns for next step
combined_anno <- combined_anno %>% tibble::column_to_rownames("Kids_First_Biospecimen_ID")
cluster_coseq <- cluster_coseq %>% tibble::column_to_rownames("gene_symbol")

# generate color list for heatmaps
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# select n distinct colors for each annotation
n_cluster_by_coseq <- cluster_coseq %>% pull(cluster_by_coseq) %>% unique() %>% length()
n_cluster_nb <- combined_anno %>% pull(cluster_assigned_nb) %>% unique() %>% length()
n_mol_subtype <- combined_anno %>% pull(molecular_subtype) %>% unique() %>% length()

# generate a list of colors for each annotation 
coseq_color <-sample(col_vector, n_cluster_by_coseq)
names(coseq_color) <- levels(cluster_coseq$cluster_by_coseq)

nb_mclust_color <- sample(col_vector, n_cluster_nb)
names(nb_mclust_color) <- levels(combined_anno$cluster_assigned_nb)

mol_subtype_color <- sample(col_vector, n_mol_subtype)
names(mol_subtype_color) <- levels(combined_anno$molecular_subtype)

anno_colors <- list(mol_subtype_color,
                    coseq_color,
                    nb_mclust_color)
names(anno_colors) <- c("molecular_subtype", "cluster_by_coseq", "cluster_assigned_nb")

# define breaks for heatmap
breaks<- seq(-3,3, length.out=92)
breaks <- c(-7,-6,-5,-4,breaks,4, 5, 6, 7)

# plot heatmap
filtered_log_cg_coding_pval_zscored_anno %>% as.matrix() %>% t() %>%
  pheatmap::pheatmap(annotation_col = combined_anno,
                     annotation_row = cluster_coseq,
                     annotation_colors = anno_colors,
                     breaks=breaks,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     cluster_rows=FALSE, 
                     cluster_cols=FALSE,
                     show_colnames = F,
                     show_rownames=F,
                     filename = file.path(output_dir, "nb_mclust_coseq_molsubtype_heatmap.pdf"), 
                     width = 10, height = 8)

# also print out a df with all samples within the cancer group and their cluster assignment
combined_anno %>% 
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>% 
  write_tsv(file.path(output_dir, "cancer_group_of_interest_nb_cluster_assigned.tsv"))
