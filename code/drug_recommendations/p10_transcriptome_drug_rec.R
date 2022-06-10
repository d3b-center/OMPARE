# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "drug_recommendations")
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "transcriptome_drug_rec.R"))

# input files
genes_up <- read.delim(file.path(patient_dir, "output", "rnaseq_analysis", 'genes_up.txt'))

# transcriptomic drug recommendations
# full output save as tsv
full_fname <- file.path(output_dir, "CHEMBL_drug-gene.tsv")
transcriptome_drug_rec_output <- transcriptome_drug_rec(diffexpr_genes = genes_up)
write.table(x = transcriptome_drug_rec_output, file = full_fname, quote = F, sep = "\t", row.names = F)

# annotate CNS penetrance from an ML-based CNS-penetrance predictor (admetSAR)
admetsar_output <- read.delim(file.path(root_dir, "data", "cmdc201800533-sup-0001-misc_information_s7.txt"))
transcriptome_drug_rec_output <- transcriptome_drug_rec_output %>%
  mutate(CNS_penetrant = ifelse(DrugBank_ID %in% admetsar_output$DrugBank.ID, "Yes", "No"))

# annotate GSNCA hub gene information to annotation
gsnca_results_dir <- file.path(patient_dir, "output", "gsnca_analysis")
gsnca_result_normal<- readr::read_tsv(file.path(gsnca_results_dir, 
                                                "cluster_samples_vs_normal_data_GSNCA_analysis.tsv"))
gsnca_result_pediatric <- readr::read_tsv(file.path(gsnca_results_dir, 
                                             "cluster_samples_vs_pediatric_data_GSNCA_analysis.tsv"))
gsnca_result_adult <- readr::read_tsv(file.path(gsnca_results_dir, 
                                               "cluster_samples_vs_adult_data_GSNCA_analysis.tsv"))

# annotate GSNCA hub gene
transcriptome_drug_rec_output <- transcriptome_drug_rec_output %>% 
  mutate(hub_GSNCA = case_when(
    grepl("Normal", Comparison) & Gene %in% gsnca_result_normal$hub_gene ~ "Yes",
    grepl("Adult", Comparison) & Gene %in% gsnca_result_adult$hub_gene ~ "Yes",
    grepl("Pediatric", Comparison) & Gene %in% gsnca_result_pediatric$hub_gene ~ "Yes",
    TRUE ~ "No"
  )
)

# subset columns for report
fname <- file.path(output_dir, "transcriptome_drug_rec.rds")
transcriptome_drug_rec_output <- transcriptome_drug_rec_output %>%
  dplyr::select(Drug, Gene, Source, Comparison, logFC, MOA, CNS_penetrant, hub_GSNCA) %>%
  unique()
saveRDS(transcriptome_drug_rec_output, file = fname)
