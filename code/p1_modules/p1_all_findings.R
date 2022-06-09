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

# source functions
source(file.path(module_dir, "utils", 'all_findings.R'))
source(file.path(module_dir, "utils", "annotate_mutations.R"))

# input maf file
patient_maf <- list.files(path = file.path(patient_dir, "simple-variants"), pattern = "consensus", full.names = T)
if(length(patient_maf) > 0){
  patient_maf <- data.table::fread(patient_maf)
  annotated_maf <- annotate_mutations(myMutData = patient_maf, myCancerGenes = cancer_genes, cancer_hotspots = cancer_hotspots)
} else {
  annotated_maf <- NULL
}

# output
fname <- file.path(output_dir, paste0("all_findings_output_", snv_caller, ".rds"))
all_findings_output <- all_findings(sample_info = sample_info,
                                    annotated_maf = annotated_maf, 
                                    filtered_fusions = filtered_fusions, 
                                    filtered_cnv = filtered_cnv, 
                                    tpm_data = tpm_data, 
                                    snv_caller = snv_caller)
saveRDS(all_findings_output, file = fname)

