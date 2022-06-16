# all findings table

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "p1_modules")
output_dir <- file.path(patient_dir, "output")
dir.create(output_dir, recursive = T, showWarnings = F)

# reference data to annotate mutations
# cancer hotspots
cancer_hotspots <- readRDS(file.path(data_dir, "cancer_hotspots.rds"))
# cancer genes
cancer_genes <- readRDS(file.path(data_dir, "cancer_gene_list.rds"))

# input data
if(file.exists(file.path(output_dir, "rnaseq_analysis", "rnaseq_analysis_output.rds"))){
  rnaseq_analysis_output <- readRDS(file.path(output_dir, "rnaseq_analysis", "rnaseq_analysis_output.rds"))
} else {
  rnaseq_analysis_output <- NULL
}

# source functions
source(file.path(module_dir, "utils", 'all_findings.R'))
source(file.path(module_dir, "utils", "annotate_mutations.R"))

# input maf file
maf_file <- list.files(path = file.path(patient_dir, "simple-variants"), pattern = "consensus", full.names = T)
if(file.exists(maf_file)){
  full_maf <- data.table::fread(maf_file, skip = 1)
  annotated_maf <- annotate_mutations(myMutData = full_maf, myCancerGenes = cancer_genes, cancer_hotspots = cancer_hotspots)
} else {
  annotated_maf <- NULL
}

# output
fname <- file.path(output_dir, "all_findings_output.rds")
all_findings_output <- all_findings(sample_info = sample_info,
                                    annotated_maf = annotated_maf, 
                                    filtered_fusions = filtered_fusions, 
                                    filtered_cnv = filtered_cnv, 
                                    rnaseq_analysis_output = rnaseq_analysis_output)
saveRDS(all_findings_output, file = fname)

