# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")
output_dir <- file.path(patient_dir, "output", "transcriptomically_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", 'pathway_analysis.R'))
source(file.path(module_dir, "utils", 'shared_pathways_plot.R'))

# load inputs
nn_table <- file.path(output_dir, "pediatric_all_nn_table.rds")
nn_table <- readRDS(nn_table)

# patient of interest
patient_of_interest <- sample_info %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  .$Kids_First_Biospecimen_ID

# pathway enrichment top correlated samples
fname <- file.path(output_dir, "pathway_analysis_pediatric.rds")
pathway_analysis_pediatric <- pathway_analysis(nn_table = nn_table, 
                                               normal_tissue = "Brain", 
                                               adult_cancer = "GBM", 
                                               pediatric_cancer = "HGAT",
                                               filtered_cnv = filtered_cnv,
                                               output_dir = output_dir, 
                                               prefix = "pediatric", 
                                               comparison = "pediatric",
                                               patient_of_interest = patient_of_interest)
# save output
saveRDS(pathway_analysis_pediatric, file = fname)

# convert to plot
shared_pathways_plot(pathway_analysis = pathway_analysis_pediatric, prefix = "pediatric", output_dir = output_dir)
