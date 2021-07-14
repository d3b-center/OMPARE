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

# Filter the histoogies file to contain only the specimens that have mutect2 results
pbta_mutect2_samples_histology <- histologies %>% 
  dplyr::filter(Kids_first_Biospecimen_ID %in% pbta_mutect2_samples)
# See what are the cohorts and experimental strategy
pbta_mutect2_samples_histology$cohort %>% unique()
pbta_mutect2_samples_histology$experimental_strategy %>% unique()

# This part pull in the bed files corresponding to the kit used in pbta
# Read in the bed files
bed_file_path = '../references'
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
bed_name_length_matched <- data.frame(bed_file_list, total_length_list)

# since the mutect file does not have PNOC008 patient, we assume PNOC is PNOC003 
pbta_mutect2_samples_histology <- pbta_mutect2_samples_histology %>%
  dplyr::mutate(bed_selected = case_when(
    experimental_strategy == "WGS" ~ "CCDS.bed",
    experimenta_strategy == "WXS" ~ "Strexome_targets_intersect_sorted_padded100.GRCh38.withCCDS.bed",
    experimenta_strategy == "Targeted Sequencing" & cohort == "PNOC" ~ "StrexomeLite_hg38_liftover_100bp_padded.bed",
    TRUE ~ "xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.withCCDS.bed"
  ))

# annotate bed_lengths to histologies - since there are only 4 bed files - do the following
pbta_mutect2_samples_histology <- pbta_mutect2_samples_histology %>%
  dplyr::mutate(bed_length = case_when(
    bed_selected == bed_name_length_matched[1,1] ~ as.numeric(bed_name_length_matched[1,2]),
    bed_selected == bed_name_length_matched[2,1] ~ as.numeric(bed_name_length_matched[2,2]),
    bed_selected == bed_name_length_matched[3,1] ~ as.numeric(bed_name_length_matched[3,2]), 
    bed_selected == bed_name_length_matched[4,1] ~ as.numeric(bed_name_length_matched[4,2])
  ))

# filter pbta
pbta_maf <- pbta_mutect2 %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
  filter(Variant_Classification %in% var_class,
         t_depth >= tumor_depth,
         vaf >= vaf_cutoff,
         t_alt_count >= var_count) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)

pbta_tmb_combined <- data.frame()
for (i in length(bed_file_list)){
  bed_name = bed_file_list[i]
  bed_length <- bed_name_length_matched %>% dplyr::filter(bed_file_list == bed_name) %>%
    dplyr::pull(total_length_list) %>% as.numeric()
  
  sample_list <- pbta_mutect2_samples_histology %>% dplyr::filter(bed_selected == bed_name) %>%
    pull(Kids_First_Biospecimen_ID) %>% unique()
  
  bed_file <- data.table::fread(file.path(paste0("../results/bed_files/tcga_not_in_pbta", bed_name)))
  colnames(bed_file)  <- c('chr', 'start', 'end')
  
  pbta_maf_matched <- pbta_maf %>% dplyr::filter(Tumor_Sample_Barcode %in% sample_list)
  
  # intersect with bed file
  subject <- with(bed_file, GRanges(chr, IRanges(start = start, end = end)))
  query <- with(pbta_maf_matched, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
  pbta_tmb_matched <- findOverlaps(query = query, subject = subject, type = "within")
  pbta_tmb_matched <- data.frame(pbta_maf_matched[queryHits(pbta_tmb),], bed_file[subjectHits(pbta_tmb),])
  
  # mutations per sample
  pbta_tmb_matched <- pbta_tmb_matched %>%
    mutate(sample_name =  Tumor_Sample_Barcode) %>%
    group_by(sample_name) %>%
    mutate(num_var = n()) %>%
    mutate(tmb = num_var*1000000/bed_length) %>%
    dplyr::select(sample_name, tmb) %>%
    unique()
  
  pbta_tmb_combined <- rbind(pbta_tmb_combined, pbta_tmb_matched)
}
    
# add short histology, bed selected and bed length
pbta_mutect2_samples_histology <- pbta_mutect2_samples_histology %>%
  dplyr::select(Kids_First_Biospecimen_ID, short_histology, bed_length) %>%
  unique()
pbta_tmb_combined <- pbta_tmb_combined %>%
  inner_join(pbta_mutect2_samples_histology, by = c('sample_name'= 'Kids_First_Biospecimen_ID')) %>%
  dplyr::select(short_histology, sample_name, tmb, bed_length)
colnames(pbta_tmb_combined) <- c("Diseasetype", "Samplename", "TMBscore", "BedLength")

# update file with new filters
write.table(pbta_tmb, file = file.path(ref_dir, 'PBTA-TMBscores_withdiseastype.txt'), quote = F, sep = "\t", row.names = F)

### plot for PNOC008-25
# # plot code from OMPARE
# tmb.calculate <- function(myTMB = tmb_bed_file, var_class = c('Missense_Mutation', 
#                                                               'Nonsense_Mutation', 
#                                                               'Frame_Shift_Del', 
#                                                               'Frame_Shift_Ins',  
#                                                               'In_Frame_Del', 
#                                                               'In_Frame_Ins'), 
#                           vaf_cutoff = 0.05, var_count = 3, tumor_depth = 25) {
#   
#   # read mutect2 for TMB profile
#   somatic.mut.pattern <- '*.maf'
#   mutFiles <- list.files(path = 'results/PNOC008-25/', pattern = somatic.mut.pattern, recursive = TRUE, full.names = T)
#   mutFiles <- grep("mutect2", mutFiles, value = TRUE)
#   if(length(mutFiles) >= 1){
#     mutFiles <- lapply(mutFiles, data.table::fread, skip = 1, stringsAsFactors = F)
#     mutData.mutect2 <- data.table::rbindlist(mutFiles)
#     mutData.mutect2 <- as.data.frame(mutData.mutect2)
#     mutData.mutect2 <- unique(mutData.mutect2)
#   }
#   
#   # apply filters using Friends of Cancer Research Project Standards 
#   myMutData <- mutData.mutect2 %>%
#     mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
#     filter(Variant_Classification %in% var_class,
#            t_depth >= tumor_depth,
#            vaf >= vaf_cutoff,
#            t_alt_count >= var_count) %>%
#     dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
#   
#   # intersect with bed file
#   subject <- with(myTMB, GRanges(chr, IRanges(start = start, end = end)))
#   query <- with(myMutData, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
#   res <- findOverlaps(query = query, subject = subject, type = "within")
#   res <- data.frame(myMutData[queryHits(res),], myTMB[subjectHits(res),])
#   
#   # return the number of variants overlapping with the bed file
#   return(nrow(res))
# }
# 
# 
# tmb_profile <- function(pedTMBScores = ped_tmb, adultTMBScores = adult_tmb, TMB, tmb_bed_file) {
#   
#   TMB <- tmb.calculate(myTMB = tmb_bed_file)/TMB
#   pedTMBScores$Type <- "Pediatric"
#   adultTMBScores$Type <- "Adult"
#   tmbScores <- rbind(pedTMBScores, adultTMBScores)
#   
#   # count median per histology
#   disease.order <- tmbScores %>%
#     group_by(Type, Diseasetype) %>%
#     summarise(median = median(TMBscore), count = n()) %>%
#     filter(count > 2) %>%
#     arrange(desc(median)) %>%
#     .$Diseasetype
#   
#   # subset to disease type with > 2 samples
#   tmbScores <- tmbScores %>%
#     filter(Diseasetype %in% disease.order)
#   tmbScores$Diseasetype <- factor(tmbScores$Diseasetype, levels = disease.order)
#   
#   # Plot it
#   print(TMB)
#   p <- ggplot(tmbScores, aes(Diseasetype, TMBscore, fill = Type)) + 
#     geom_boxplot() + theme_bw() + scale_y_log10(breaks = c(.25, 1, 10, 100, 500)) + 
#     scale_fill_manual(values = c("blue", "red")) + 
#     xlab("Disease") + ylab("Mutations per MB") + 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#     geom_hline(yintercept = TMB, linetype = 2, color = 'gray30') +
#     annotate("text", x = 35, y = max(tmbScores$TMBscore) - 50, 
#              label = "- - - Patient TMB", size = 4, 
#              fontface = 'italic', color = "gray30")
#   return(p)
# }
# 
# # read updated files
# ped_tmb <- data.table::fread(file.path(ref_dir, 'pbta-TMBscores_withdiseastype.txt'))
# adult_tmb <- data.table::fread(file.path(ref_dir, 'TCGA_diseasetypes_and_samples_TMBscores.txt'))
# tmb_bed_file <- data.table::fread(file.path(ref_dir, 'xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed'))
# colnames(tmb_bed_file)  <- c('chr', 'start', 'end')
# tmb = 77.46
# tmb_profile(pedTMBScores = ped_tmb, 
#             adultTMBScores = adult_tmb, 
#             TMB = tmb, 
#             tmb_bed_file = tmb_bed_file)

