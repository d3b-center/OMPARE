library(tidyverse)
library(GenomicRanges)

# define scratch directory
scratch_dir = '/tmp/'

# set root directory and other directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
results_dir <- file.path(root_dir, 'tmb_normalization', 'results')
ref_dir <- file.path(root_dir, 'tmb_normalization', 'references')

# filters
var_class = c('Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',  'In_Frame_Del', 'In_Frame_Ins')
vaf_cutoff = 0.05
var_count = 3
tumor_depth = 25

# function to merge gdc downloaded mutect2 maf files into 1 data frame and add the disease_type 
merge_files <- function(nm){
  disease_type <- gsub(".*TCGA[.]", "", nm)
  disease_type <- gsub("[.].*", "", disease_type)
  x <- data.table::fread(nm)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    x$disease_type <- disease_type
    return(x)
  }
}

# Read in the list of maf files and get the tumor type 
tcga_mutect2 <- list.files(path = file.path(scratch_dir, 'maf_files'), pattern = '.maf.gz', recursive = TRUE, full.names = T)
tcga_mutect2 <- lapply(tcga_mutect2, FUN = function(x) merge_files(x))
tcga_mutect2 <- data.table::rbindlist(tcga_mutect2, fill = T)

# Read in the manifest for bed files
tcga_bed_manifest <- read.delim(file.path(results_dir, "tcga_not_in_pbta_bam_manifest_with_length.tsv")) %>%
  dplyr::rename(Tumor_Sample_Barcode = associated_entities.entity_submitter_id) %>%
  dplyr::select(Tumor_Sample_Barcode, bed_selected, bed_length)

tcga_bed_selection <- read.delim(file.path(results_dir, "tcga_not_in_pbta_unique_bed_with_length.tsv"), header = T, sep = "\t", stringsAsFactor = F)

# find out the unique bed files
tcga_bed_list <- tcga_bed_selection$bed_selected %>% unique() 

# calculate the maf 
tcga_maf_combined <- tcga_mutect2 %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
  filter(Variant_Classification %in% var_class,
         t_depth >= tumor_depth,
         vaf >= vaf_cutoff,
         t_alt_count >= var_count) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position, disease_type, Tumor_Sample_Barcode)

tcga_maf_combined <- tcga_maf_combined %>% left_join(tcga_bed_manifest)

# do the intersection of maf files and the bed files using corresponding bed files
tcga_tmb_list <- lapply(tcga_bed_list, function(x){
    # read in the correct bed file
    bed_file <- data.table::fread(file.path(results_dir, 'bed_files', 'tcga_not_in_pbta', x))
    bed_file <- bed_file[,1:3]
    colnames(bed_file)  <- c('chr', 'start', 'end')
    bed_length <- tcga_bed_selection %>% 
      dplyr::filter(bed_selected == x) %>%
      pull(bed_length) %>% as.numeric()
    
    # find the maf results of the sample 
    tcga_maf <- tcga_maf_combined %>% dplyr::filter(bed_selected == x)
  
    # intersect with bed file
    subject <- with(bed_file, GRanges(chr, IRanges(start = start, end = end)))
    query <- with(tcga_maf, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
    tcga_tmb <- findOverlaps(query = query, subject = subject, type = "within")
    tcga_tmb <- data.frame(tcga_maf[queryHits(tcga_tmb),], bed_file[subjectHits(tcga_tmb),])
  
    # mutations per sample
    tcga_tmb <- tcga_tmb %>%
      mutate(sample_name =  Tumor_Sample_Barcode) %>%
      group_by(sample_name) %>%
      mutate(num_var = n()) %>%
      mutate(tmb = num_var*1000000/bed_length) %>%
      dplyr::select(disease_type, sample_name, tmb) %>%
      unique() %>% mutate(BedLength = bed_length)
    colnames(tcga_tmb) <- c("Diseasetype", "Samplename", "TMBscore", "BedLength")
  return(tcga_tmb)
  }
)
tcga_tmb_combined <- do.call(rbind,tcga_tmb_list)

# update file with new filters and save to data directory
write.table(tcga_tmb_combined, file = file.path(data_dir, 'tmb', 'TCGA_not_in_pbta_diseasetypes_and_samples_TMBscores.txt'), quote = F, sep = "\t", row.names = F)


##############################################################################
# The following code generate the output for tcga samples that are included
# in the pbta project
##############################################################################

# Read in the list of maf files and get the tumor type 
tcga_mutect2 <- data.table::fread(file.path(ref_dir, 'pbta-tcga-snv-mutect2.vep.maf.gz'))

# Read in the manifest for bed files
tcga_bed_manifest <- read.delim(file.path(results_dir, "tcga_in_pbta_bam_manifest_with_length.tsv")) %>%
  dplyr::select(Tumor_Sample_Barcode, bed_selected, bed_length, Primary_diagnosis)
tcga_bed_list <- tcga_bed_manifest$bed_selected %>% unique()

# calculate the maf 
tcga_maf_combined <- tcga_mutect2 %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
  filter(Variant_Classification %in% var_class,
         t_depth >= tumor_depth,
         vaf >= vaf_cutoff,
         t_alt_count >= var_count) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)

# add the bed file and bed length to the file 
tcga_maf_combined <- tcga_maf_combined %>% left_join(tcga_bed_manifest)

# do the intersection of maf files and the bed files using corresponding bed files
tcga_tmb_list <- lapply(tcga_bed_list, function(x){
  # find the name of the bed file 
  bed_file <- data.table::fread(file.path(results_dir, 'bed_files', 'tcga_in_pbta', x))
  bed_file <- bed_file[,1:3]
  colnames(bed_file)  <- c('chr', 'start', 'end')
  bed_length <- tcga_bed_manifest %>% dplyr::filter(bed_selected == x) %>%
    pull(bed_length) %>% unique() %>% as.numeric()
  
  # find the maf results of the sample 
  tcga_maf <- tcga_maf_combined %>% dplyr::filter(bed_selected == x)
  
  # intersect with bed file
  subject <- with(bed_file, GRanges(chr, IRanges(start = start, end = end)))
  query <- with(tcga_maf, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
  tcga_tmb <- findOverlaps(query = query, subject = subject, type = "within")
  tcga_tmb <- data.frame(tcga_maf[queryHits(tcga_tmb),], bed_file[subjectHits(tcga_tmb),])
  
  # mutations per sample
  tcga_tmb <- tcga_tmb %>%
    mutate(sample_name =  Tumor_Sample_Barcode) %>%
    group_by(sample_name) %>%
    mutate(num_var = n()) %>%
    mutate(tmb = num_var*1000000/bed_length) %>%
    dplyr::select(Primary_diagnosis, sample_name, tmb) %>%
    unique() %>% mutate(BedLength = bed_length)
  colnames(tcga_tmb) <- c("Diseasetype", "Samplename", "TMBscore", "BedLength")
  return(tcga_tmb)
}
)
tcga_tmb_combined <- do.call(rbind,tcga_tmb_list)

# update file with new filters and save to data directory
write.table(tcga_tmb_combined, file = file.path(data_dir, 'tmb', 'TCGA_in_pbta_diseasetypes_and_samples_TMBscores.txt'), quote = F, sep = "\t", row.names = F)
