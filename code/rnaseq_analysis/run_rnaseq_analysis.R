# Author: Komal S. Rathi
# Function: differential genes + pathways for patient of interest sample vs GTEx normal tissues

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "rnaseq_analysis")
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions for RNA-seq diffexpr & pathway analysis
source(file.path(module_dir, "utils", "rnaseq_analysis_edgeR.R"))

# gene set (c2 reactome)
gene_set <- getGmt(file.path(data_dir, 'msigdb', 'c2.cp.reactome.v6.0.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
gene_set <- geneIds(gene_set)
gene_set_ts <- stack(gene_set)

# cancer genes
cancer_genes <- readRDS(file.path(root_dir, "data", "cancer_gene_list.rds"))

# gtex tissue counts
gtex_counts <- list.files(path = file.path("data", "normal_data"), pattern = "*counts.rds", full.names = T)
gtex_counts <- readRDS(gtex_counts)

# format input tpm data
tpm_data_melted <- tpm_data %>% 
  rownames_to_column("gene_symbol") %>%
  gather(sample, tpm, -c('gene_symbol'))

# format input count data
count_data_melted <- count_data %>% 
  rownames_to_column("gene_symbol") %>%
  gather(sample, expected_count, -c('gene_symbol'))

# name of comparison
patient_of_interest <- sample_info %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  .$Kids_First_Biospecimen_ID
comparison <- paste0(patient_of_interest, '_vs_GTExBrain')

# module directory
output_dir <- file.path(patient_dir, "output", "rnaseq_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# output file
fname <- file.path(output_dir, "rnaseq_analysis_output.rds")

# run single sample level edgeR analysis
rnaseq_analysis_output <- run_rnaseq_analysis_edger(exp.data.tpm = tpm_data_melted, # tpm
                                                    exp.data.counts = count_data_melted, # counts
                                                    refData.counts = gtex_counts, # gtex counts
                                                    gene_set = gene_set, # reactome
                                                    comparison = comparison, # comparison name
                                                    single_sample = TRUE,    # single sample analysis
                                                    sample_name = patient_of_interest, # pass sample name only for single sample analysis 
                                                    cancer_genes = cancer_genes)   

# save output
saveRDS(rnaseq_analysis_output, file = fname)
