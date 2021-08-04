BiocManager::install("GenomicDataCommons")

suppressPackageStartupMessages({
  library(TCGAutils)
  library(jsonlite)
  library(curl)
  library(downloader)
  library(GenomicDataCommons) # use GenomicDataCommons package to generate manifest file
  library(magrittr)
})

## View available data through GenomicDataCommons
available_values('files','cases.project.project_id')
available_values('files','experimental_strategy')
available_values('files','data_format')

# read in the previous TCGA tumor score file to find out the TCGA projects used:
tcga_previous <- read.delim("../references/TCGA_diseasetypes_and_samples_TMBscores.txt", header = T, sep = "\t") %>%
  dplyr::mutate(tcga_project_id = paste0("TCGA-", Diseasetype))
tcga_project_list <- tcga_previous %>% dplyr::pull(tcga_project_id) %>% unique()

## Define fQuery_maf for subsequent query from GDC portal 
fQuery_maf = files()
default_fields(fQuery_maf)

## Query GDC for MAF files using experimental strategy, project id and data format
fQuery.maf.files = GenomicDataCommons::filter(fQuery_maf,~ cases.project.project_id == tcga_project_list &
                                              experimental_strategy %in% c("WXS", "WGS", "Targeted Sequencing") &
                                              data_format == "MAF") %>%
  GenomicDataCommons::manifest()

## Subset it to mutect file containing somatic maf files 
fQuery_maf_list <- fQuery.maf.files %>% dplyr::filter(grepl("mutect", filename)) %>% 
  dplyr::filter(grepl("somatic.maf.gz", filename)) %>%
  pull(id)

# Find the UUID to barcode for all the specimens in the file 
maf_manifest_list <- lapply(fQuery_maf_list, function(x){
  fQuery.each.maf <- fQuery.maf.files %>% dplyr::filter(id == x) 
  maf_manifest_each <- cbind(fQuery.each.maf, UUIDtoBarcode(fQuery.each.maf$id, from_type = "file_id"))
  return(maf_manifest_each)
})
maf_manifest <- do.call(rbind, maf_manifest_list)

# Decipher the tumor sample barcode 
maf_manifest <- cbind(maf_manifest, TCGAutils::TCGAbiospec(maf_manifest$associated_entities.entity_submitter_id))
# Get the disease type from the name of the maf file
maf_manifest <- maf_manifest %>% dplyr::mutate(disease_type = gsub(".*TCGA[.]", "", filename)) %>%
  dplyr::mutate(disease_type = gsub("[.].*", "",disease_type))
  
## Define fQuery_bam for subsequent query from GDC portal 
fQuery_bam = files()
default_fields(fQuery_bam)

## Query GDC for BAM files using experimental strategy, project id and data format
fQuery.bam.files = GenomicDataCommons::filter(fQuery_bam,~ cases.project.project_id == tcga_project_list &
                                              experimental_strategy %in% c("WXS", "WGS", "Targeted Sequencing") &
                                              data_format == "BAM") %>%
  GenomicDataCommons::manifest()

# Find the list of files 
fQuery_bam_list<-fQuery.bam.files %>% dplyr::pull(id)

# Find the UUID to barcode for all the specimens in the file 
bam_manifest_list <- lapply(fQuery_bam_list, function(x){
  fQuery.each.bam <- fQuery.bam.files %>% dplyr::filter(id == x) 
  bam_manifest_each <- cbind(fQuery.each.bam, UUIDtoBarcode(fQuery.each.bam$id, from_type = "file_id"))
  return(bam_manifest_each)
})
bam_manifest <- do.call(rbind, bam_manifest_list)

# check whether all maf files have bam files
all(maf_manifest$associated_entities.entity_submitter_id %in% bam_manifest$associated_entities.entity_submitter_id)
# check whether all bam files have maf files
all(bam_manifest$associated_entities.entity_submitter_id %in% maf_manifest$associated_entities.entity_submitter_id)

# Not all bam files have maf files and for our purposes we want to correct TMB so only interested 
# in files with both bam (bed) and mutect2 (maf) files

bam_manifest <- bam_manifest %>% 
  dplyr::filter(bam_manifest$associated_entities.entity_submitter_id %in% maf_manifest$associated_entities.entity_submitter_id)
bam_manifest_list<-bam_manifest %>% pull(id) %>% unique()

## curl information about target capture kits
res <- lapply(bam_manifest_list, function(x) {
  con = curl::curl(paste0("https://api.gdc.cancer.gov/files/", 
                          x, "?pretty=true&fields=analysis.metadata.read_groups.target_capture_kit_target_region,analysis.metadata.read_groups.target_capture_kit_name"))
  x = jsonlite::fromJSON(con)
  return(x)
})

# Grabbing the target_capture_kit_name, target_capture_kit_vendor, target_capture_kit_catalog_number and target_capture_kit_target_region
# from the data.analysis.metadata.read_groups return of the GDC curl
request_list <- lapply(res, function(x) unique(x$data$analysis$metadata$read_groups))
request_list <- do.call(rbind, request_list)
bam_manifest <- cbind(bam_manifest, request_list)

# there are also files that we do not know which bed files are used so we will just go ahead and delete those
bam_manifest <- bam_manifest %>% dplyr::filter(!is.na(target_capture_kit_target_region))
maf_manifest <- maf_manifest %>% 
  dplyr::filter(maf_manifest$associated_entities.entity_submitter_id %in% bam_manifest$associated_entities.entity_submitter_id)

# read out the unique bed files that are used
unique_bed_files <- bam_manifest %>% dplyr::select(target_capture_kit_name,target_capture_kit_target_region) %>% 
  unique() %>% data.frame()

# write out all the tsv files for the next steps
readr::write_tsv(bam_manifest, "../results/tcga_not_in_pbta_bam_manifest.tsv")
readr::write_tsv(maf_manifest, "../results/tcga_not_in_pbta_maf_manifest.tsv")
readr::write_tsv(unique_bed_files, "../results/tcga_not_in_pbta_unique_bed_file.tsv")

