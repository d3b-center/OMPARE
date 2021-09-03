# all findings table

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "p1_modules")

# source functions
source(file.path(module_dir, "utils", 'all_findings.R'))

# output
all_findings_output <- all_findings(annotated_mutations = mutDataAnnot, 
                                    filtered_fusions = fusData, 
                                    filtered_cnvs = cnvDataFilt, 
                                    expr_data = expData, 
                                    snv_caller = snv_caller)

