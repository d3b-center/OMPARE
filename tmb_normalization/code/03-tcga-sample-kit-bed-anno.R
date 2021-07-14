# This script uese all the downloaded bed files, calculate the length the bed files
# And map the bed length back to the bam_manifest file

# Read in the bed files
bed_file_path = '../results/bed_files/tcga_not_in_pbta'
bed_file_list <- list.files(path = bed_file_path, pattern = '.bed', recursive = TRUE, full.names = T)
bed_files <- lapply(bed_file_list, read.delim, header=F)
# Now convert the bed file list to ones only contain the bed file name
bed_file_list <- gsub(paste0(bed_file_path, "/"), "", bed_file_list)

# Read in the manifest file and get a column that only contains the bed file (name Gh38 corrected)
bam_manifest <- read.delim("../results/tcga_not_in_pbta_bam_manifest.tsv")

# for some of the kits, multiple target kits are listed since GDC is not sure which one is used
# some of those kits do not have associated bed files - so we just use the files that I generated 
# before and select one fothe the target kits that actually have a bed file 
bed_selection <- read.delim("../results/tcga_not_in_pbta_bed_selected.tsv", header = T, row.names=NULL, stringsAsFactors = F)
bam_manifest <- bam_manifest %>% dplyr::left_join(bed_selection)

# Define a vector that will hold the total lengths for all the bed files in the bed file list
total_length_list <- c()
# Calculate bed region length for the bed file
for (j in 1:length(bed_file_list)){
  total_length<-0
  bed_file <- bed_files[[j]]
  for (i in 1:nrow(bed_file)) {
    distance <- as.numeric(bed_file[i,3]) - as.numeric(bed_file[i,2])
    total_length = total_length + distance
    total_length_list[j] = total_length
  }
}

# Generate a dataframe containing the bed file name and the bed length
bed_name_length_matched <- data.frame(bed_file_list, total_length_list) %>%
  dplyr::rename(bed_selected = bed_file_list, bed_length = total_length_list)
# match the bed length to the manifest
bam_manifest <- bam_manifest %>% dplyr::left_join(bed_name_length_matched)

# Write out the updated bam_manifest file
readr::write_tsv(bam_manifest, "../results/tcga_not_in_pbta_bam_manifest_with_length.tsv")
readr::write_tsv(bed_name_length_matched, "../results/tcga_not_in_pbta_unique_bed_with_length.tsv")

##############################################################################
## The following code calculate the bed length for tcga in pbta 
##############################################################################

# Read in the bed files
bed_file_path = '../results/bed_files/tcga_in_pbta'
bed_file_list <- list.files(path = bed_file_path, pattern = '.Gh38.bed', recursive = TRUE, full.names = T)
bed_files <- lapply(bed_file_list, read.delim, header=F)
# Now convert the bed file list to ones only contain the bed file name
bed_file_list <- gsub(paste0(bed_file_path, "/"), "", bed_file_list)

# Read in the manifest file and get a column that only contains the bed file (name Gh38 corrected)
bam_manifest <- read.delim("../references/pbta-tcga-manifest.tsv")

# for the next steps, since we need to use the tumor sample barcode to match Tumor_Sample_Barcode 
# in the mutect2 file, we need to generate another column containing the tumor sample barcode
bam_manifest <- bam_manifest %>% dplyr::mutate(Tumor_Sample_Barcode = gsub(".*[.]TCGA", "TCGA", Tumor_BAM)) %>%
  dplyr::mutate(Tumor_Sample_Barcode = gsub("[.].*", "", Tumor_Sample_Barcode))

# one thing is that previously, Yuankun had been taking the intersect of all the listed bed files
# but to be consistent with this TCGA workflow, I only choose the first listed bed file 
bam_manifest <- bam_manifest %>% dplyr::mutate(bed_selected = dplyr::case_when(
  BED_In_Use == "intersected_whole_exome_agilent_designed_120_AND_tcga_6k_genes.Gh38.bed" ~ "tcga_6k_genes.targetIntervals.Gh38.bed", 
  BED_In_Use == "intersected_whole_exome_agilent_plus_tcga_6k_AND_tcga_6k_genes.Gh38.bed" ~ "whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed", 
  BED_In_Use == "whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.Gh38.bed" ~ "whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.Gh38.bed",
  BED_In_Use == "whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed" ~ "whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed",
  TRUE ~ "whole_exome_agilent_designed_120.targetIntervals.Gh38.bed"
  )
)

# for tcga in pbta, the bed file in use is already annotated to BED_In_Use field
# Define a vector that will hold the total lengths for all the bed files in the bed file list
total_length_list <- c()
# Calculate bed region length for the bed file
for (j in 1:length(bed_file_list)){
  total_length<-0
  bed_file <- bed_files[[j]]
  for (i in 1:nrow(bed_file)) {
    distance <- as.numeric(bed_file[i,3]) - as.numeric(bed_file[i,2])
    total_length = total_length + distance
    total_length_list[j] = total_length
  }
}

# Generate a dataframe containing the bed file name and the bed length
bed_name_length_matched <- data.frame(bed_file_list, total_length_list) %>%
  dplyr::rename(bed_selected = bed_file_list, bed_length = total_length_list)

# Add total length of bed file to the manifest file 
bam_manifest <- bam_manifest %>% dplyr::left_join(bed_name_length_matched)

# Write out the updated bam_manifest file
readr::write_tsv(bam_manifest, "../results/tcga_in_pbta_bam_manifest_with_length.tsv")
readr::write_tsv(bed_name_length_matched, "../results/tcga_in_pbta_unique_bed_with_length.tsv")
