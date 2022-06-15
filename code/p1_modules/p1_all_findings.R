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
maf_file <- list.files(path = file.path(patient_dir, "simple-variants"), pattern = snv_caller, full.names = T)
if(is.na(snv_caller)){
  maf_file <- NULL
  annotated_maf <- NULL
} else if(snv_caller == "all"){ 
  maf_file <- grep('consensus', maf_file, invert = TRUE, value = TRUE)
} else {
  maf_file <- grep(snv_caller, maf_file, value = TRUE)
}

if(length(maf_file) >= 1){
  maf_file <- lapply(maf_file, data.table::fread, skip = 1, stringsAsFactors = F)
  full_maf <- data.table::rbindlist(maf_file, fill = TRUE)
  full_maf <- as.data.frame(full_maf)
  full_maf <- unique(full_maf)
  annotated_maf <- annotate_mutations(myMutData = full_maf, myCancerGenes = cancer_genes, cancer_hotspots = cancer_hotspots)
} 

# output
fname <- file.path(output_dir, paste0("all_findings_output_", snv_caller, ".rds"))
all_findings_output <- all_findings(sample_info = sample_info,
                                    annotated_maf = annotated_maf, 
                                    filtered_fusions = filtered_fusions, 
                                    filtered_cnv = filtered_cnv, 
                                    rnaseq_analysis_output = rnaseq_analysis_output, 
                                    snv_caller = snv_caller)
saveRDS(all_findings_output, file = fname)

