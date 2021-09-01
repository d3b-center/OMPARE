# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# sample of interest
pnoc008_sample_of_interest <- sampleInfo$subjectID

# references
patient_dir <- file.path(root_dir, "results", pnoc008_sample_of_interest)
chembldb_path <- file.path(ref_dir, "chembl", "chembl_29_sqlite", "chembl_29.db")

# pnoc008 clinical file
pnoc008_clinical <- file.path(ref_dir, "pnoc008", "pnoc008_clinical.rds")

# inputs from patient's cemitools directory
cemitools_dir <- file.path(patient_dir, "CEMiTools")
interaction <- file.path(cemitools_dir, "interactions.tsv")
enrichment_nes <- file.path(cemitools_dir, "enrichment_nes.tsv")
cluster <- file.path(cemitools_dir, "clustered_samples.rds")

# inputs from patient's output directory (lincs connectivity analysis output)
gtex_qSig <- file.path(patient_dir, "output", "GTExBrain_qSig_output.txt")
pbta_qSig <- file.path(patient_dir, "output", "PBTA_ALL_qSig_output.txt")
pbta_hgg_qSig <- file.path(patient_dir, "output", "PBTA_HGG_qSig_output.txt")

# output files that will be saved under OMPARE/results/PNOC008-XX/output
# subnetwork file and drug mapped subnetwork qSig file
subnetwork <- file.path(patient_dir, "output", "subnetwork_genes.tsv")
subnetwork_mapped <- file.path(patient_dir, "output", "subnetwork_gene_drug_map.tsv")
gtex_subnet_qSig_mapped <- file.path(patient_dir, "output", "gtex_qSig_subnetwork_drug_gene_map.tsv")
pbta_subnet_qSig_mapped <- file.path(patient_dir, "output", "pbta_qSig_subnetwork_drug_gene_map.tsv")
pbta_hgg_subnet_qSig_mapped <- file.path(patient_dir, "output", "pbta_hgg_qSig_subnetwork_drug_gene_map.tsv")

# synergy score for all comparisons
output_gtex <- file.path(patient_dir, "output", "gtex_qSig_synergy_score.tsv")
output_pbta <- file.path(patient_dir, "output", "pbta_qSig_synergy_score.tsv")
output_pbta_hgg <- file.path(patient_dir, "output", "pbta_hgg_qSig_synergy_score.tsv")
output_combined <- file.path(patient_dir, "output", "combined_qSig_synergy_score.tsv")

# module path
drug_synergy_module <- file.path(root_dir, "code", "drug_synergy")

# run script to obtain drugs that are associated with all the genes in the subnetwork
subnetwork_qSig_gene_drug_map <- file.path(drug_synergy_module, "01-subnetwork_qSig_gene_drug_map.R")
cmd <- paste('Rscript', subnetwork_qSig_gene_drug_map,
       '--interaction', interaction,
       '--enrichment_nes', enrichment_nes,
       '--cluster', cluster,
       '--clinical', pnoc008_clinical,
       '--chemblDb_path', chembldb_path,
       '--sample_of_interest', pnoc008_sample_of_interest,
       '--gtex_qSig', gtex_qSig,
       '--pbta_qSig', pbta_qSig,
       '--pbta_hgg_qSig', pbta_hgg_qSig,
       '--subnetwork', subnetwork, 
       '--subnetwork_mapped', subnetwork_mapped,
       '--gtex_mapped', gtex_subnet_qSig_mapped,
       '--pbta_mapped', pbta_subnet_qSig_mapped,
       '--pbta_hgg_mapped', pbta_hgg_subnet_qSig_mapped)
system(cmd)

# run script to obtain drugs that are both in qSig and subnetwork
drug_synergy_score_calc <- file.path(drug_synergy_module, "02-drug_synergy_score_calc.R")
cmd <- paste('Rscript', drug_synergy_score_calc,
              '--subnetwork', subnetwork, 
              '--subnetwork_mapped', subnetwork_mapped,
              '--gtex_mapped', gtex_subnet_qSig_mapped,
              '--pbta_mapped', pbta_subnet_qSig_mapped,
              '--pbta_hgg_mapped', pbta_hgg_subnet_qSig_mapped,
              '--output_gtex', output_gtex,
              '--output_pbta', output_pbta,
              '--output_pbta_hgg', output_pbta_hgg,
              '--output_combined', output_combined)
system(cmd)
