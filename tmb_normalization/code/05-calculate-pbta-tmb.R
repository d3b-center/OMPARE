library(tidyverse)
library(GenomicRanges)

# filters
var_class = c('Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',  'In_Frame_Del', 'In_Frame_Ins')
vaf_cutoff = 0.05
var_count = 3
tumor_depth = 25

# read in the PBTA mutect2 file 
pbta_mutect2 <- data.table::fread('../references/pbta-snv-mutect2.vep.maf.gz')

# read in the histologies file
histologies <- read.delim("../references/pbta-histologies.tsv", header = T, sep = "\t", stringsAsFactors = F)

# Take a look at all the samples that are in there 
pbta_mutect2_samples <- pbta_mutect2$Tumor_Sample_Barcode %>% unique()

# Filter the histologies file to contain only the specimens that have mutect2 results
pbta_mutect2_samples_histology <- histologies %>% 
  dplyr::filter(Kids_First_Biospecimen_ID %in% pbta_mutect2_samples)
# See what are the cohorts and experimental strategy
pbta_mutect2_samples_histology$cohort %>% unique()
pbta_mutect2_samples_histology$experimental_strategy %>% unique()

# This part pull in the bed files corresponding to the kit used in pbta
# Read in the bed files
bed_file_path = '../references/pbta_bed_files'
bed_file_list <- list.files(path = bed_file_path, pattern = '.bed', recursive = TRUE, full.names = T)
bed_files <- lapply(bed_file_list, read.delim, header=F)
# Now convert the bed file list to ones only contain the bed file name
bed_file_list <- gsub(paste0(bed_file_path, "/"), "", bed_file_list)

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

# since the mutect file does not have PNOC008 patient, we assume PNOC is PNOC003 
# for some reason the `StrexomeLite_hg38_liftover_100bp_padded.bed` is no longer available - hence used
# another one to substitute

pbta_mutect2_samples_histology <- pbta_mutect2_samples_histology %>%
  dplyr::mutate(bed_selected = case_when(
    experimental_strategy == "WGS" ~ "CCDS.bed",
    experimental_strategy == "WXS" ~ "Strexome_targets_intersect_sorted_padded100.GRCh38.withCCDS.bed",
    experimental_strategy == "Targeted Sequencing" & cohort == "PNOC" ~ "Strexome_targets_intersect_sorted_padded100.GRCh38.withCCDS.bed",
    TRUE ~ "xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.withCCDS.bed"
  )) 

# annotate bed_lengths to histologies 
pbta_mutect2_samples_histology <- pbta_mutect2_samples_histology %>%
  dplyr::left_join(bed_name_length_matched) %>% 
  dplyr::mutate(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Tumor_Sample_Barcode, bed_selected, bed_length, short_histology, experimental_strategy, cohort)

# filter pbta maf file 
pbta_maf <- pbta_mutect2 %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
  filter(Variant_Classification %in% var_class,
         t_depth >= tumor_depth,
         vaf >= vaf_cutoff,
         t_alt_count >= var_count) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)

# add bed file length and bed file selected to pbta_maf
pbta_maf <- pbta_maf %>% left_join(pbta_mutect2_samples_histology)

# calculate maf results based on matching bed files
pbta_maf_list<-lapply(bed_file_list, function(x){
  bed_length <- bed_name_length_matched %>% dplyr::filter(bed_selected == x) %>%
    dplyr::pull(bed_length) %>% as.numeric()
  
  # read in the bed files for the next step overlaping
  bed_file <- data.table::fread(file.path(paste0("../references/pbta_bed_files/", x)))
  bed_file <- bed_file[,1:3]
  colnames(bed_file)  <- c('chr', 'start', 'end')
  
  # filter to only the samples that use the same bed file
  pbta_maf_matched <- pbta_maf %>% dplyr::filter(bed_selected == x)
  
  # intersect with bed file
  subject <- with(bed_file, GRanges(chr, IRanges(start = start, end = end)))
  query <- with(pbta_maf_matched, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
  pbta_tmb_matched <- findOverlaps(query = query, subject = subject, type = "within")
  pbta_tmb_matched <- data.frame(pbta_maf_matched[queryHits(pbta_tmb_matched),], bed_file[subjectHits(pbta_tmb_matched),])
  
  # mutations per sample
  pbta_tmb_matched <- pbta_tmb_matched %>%
    mutate(sample_name =  Tumor_Sample_Barcode) %>%
    group_by(sample_name) %>%
    mutate(num_var = n()) %>%
    mutate(tmb = num_var*1000000/bed_length) %>%
    dplyr::select(sample_name, short_histology, cohort, experimental_strategy, tmb, bed_length) %>%
    unique()
})

pbta_tmb_combined <- do.call(rbind,pbta_maf_list)
colnames(pbta_tmb_combined) <- c("Samplename", "Diseasetype", "Cohort", "ExpStrategy", "TMBscore", "BedLength")

# update file with new filters
write.table(pbta_tmb_combined, file = '../results/PBTA-TMBscores_withdiseastype.txt', quote = F, sep = "\t", row.names = F)

