# Author: Komal S. Rathi
# Function to create germline marker list

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

pharmacogenomics.genes <- read.delim(file.path(ref_dir, "Pharmacogenomics_Genes.list"), stringsAsFactors = F, header = F)
pharmacogenomics.genes <- data.frame("Gene" = pharmacogenomics.genes$V1, "Class" = "Pharmacogenomics")
chop.panel.genes <- read.delim(file.path(ref_dir, "CHOP_Additional_Cancer_Genes.list"), stringsAsFactors = F, header = F)
chop.panel.genes <-  data.frame("Gene" = chop.panel.genes$V1, "Class" = "CHOP Panel")
acmg.genes <- read.delim(file.path(ref_dir, "ACMG_22_Cancer_Genes.list"), stringsAsFactors = F, header = F)
acmg.genes <-  data.frame("Gene" = acmg.genes$V1, "Class" = "ACMG")
germlineMarkers <- rbind(pharmacogenomics.genes,  acmg.genes, chop.panel.genes)

# save output
saveRDS(germlineMarkers, file.path(ref_dir, "germline_markers_list.rds"))
