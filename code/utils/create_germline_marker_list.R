# Author: Komal S. Rathi
# Function to create germline marker list

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

pharmacogenomics.genes <- read.delim(file.path(data_dir, "input_files", "Pharmacogenomics_Genes.list"), stringsAsFactors = F, header = F)
pharmacogenomics.genes <- data.frame("Gene" = pharmacogenomics.genes$V1, "Class" = "Pharmacogenomics")
chop.panel.genes <- read.delim(file.path(data_dir, "input_files", "CHOP_Additional_Cancer_Genes.list"), stringsAsFactors = F, header = F)
chop.panel.genes <-  data.frame("Gene" = chop.panel.genes$V1, "Class" = "CHOP Panel")
acmg.genes <- read.delim(file.path(data_dir, "input_files" ,"ACMG_22_Cancer_Genes.list"), stringsAsFactors = F, header = F)
acmg.genes <-  data.frame("Gene" = acmg.genes$V1, "Class" = "ACMG")
germlineMarkers <- rbind(pharmacogenomics.genes,  acmg.genes, chop.panel.genes)

# save output
saveRDS(germlineMarkers, file.path(data_dir, "germline_markers_list.rds"))
