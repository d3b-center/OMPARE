library(tidyverse)
library(GenomicRanges)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
open_pbta_analysis <- '~/Projects/OpenPBTA-analysis/data/'

# bed file
tmb_bed_file <- data.table::fread(file.path(ref_dir, 'xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed'))
colnames(tmb_bed_file)  <- c('chr', 'start', 'end')

# filters
var_class = c('Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',  'In_Frame_Del', 'In_Frame_Ins')
vaf_cutoff = 0.05
var_count = 3
tumor_depth = 25

# function to merge gdc downloaded mutect2 maf files into 1 data frame
merge_files <- function(nm){
  tcga_subtype <- gsub(".*TCGA[.]", "", nm)
  tcga_subtype <- gsub("[.].*", "", tcga_subtype)
  x <- data.table::fread(nm)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    x$tcga_subtype <- tcga_subtype
    return(x)
  }
}

tcga_mutect2 <- list.files(path = '~/Downloads/gdc_download_20210208_190643.670297/', pattern = '.maf.gz', recursive = TRUE, full.names = T)
tcga_mutect2 <- lapply(tcga_mutect2, FUN = function(x) merge_files(x))
tcga_mutect2 <- data.table::rbindlist(tcga_mutect2, fill = T)

# filter tcga
tcga_maf <- tcga_mutect2 %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
  filter(Variant_Classification %in% var_class,
         t_depth >= tumor_depth,
         vaf >= vaf_cutoff,
         t_alt_count >= var_count) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position, tcga_subtype)

# intersect with bed file
subject <- with(tmb_bed_file, GRanges(chr, IRanges(start = start, end = end)))
query <- with(tcga_maf, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
tcga_tmb <- findOverlaps(query = query, subject = subject, type = "within")
tcga_tmb <- data.frame(tcga_maf[queryHits(tcga_tmb),], tmb_bed_file[subjectHits(tcga_tmb),])

# mutations per sample
tcga_tmb <- tcga_tmb %>%
  mutate(sample_name =  Tumor_Sample_Barcode) %>%
  group_by(sample_name) %>%
  mutate(num_var = n()) %>%
  mutate(tmb = num_var/77.46) %>%
  dplyr::select(tcga_subtype, sample_name, tmb) %>%
  unique()
colnames(tcga_tmb) <- c("Diseasetype", "Samplename", "TMBscore")

# update file with new filters
write.table(tcga_tmb, file = file.path(ref_dir, 'TCGA_diseasetypes_and_samples_TMBscores.txt'), quote = F, sep = "\t", row.names = F)

# PBTA
pbta_mutect2 <- data.table::fread(file.path(open_pbta_analysis, 'pbta-snv-mutect2.vep.maf.gz'))

# filter pbta
pbta_maf <- pbta_mutect2 %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
  filter(Variant_Classification %in% var_class,
         t_depth >= tumor_depth,
         vaf >= vaf_cutoff,
         t_alt_count >= var_count) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)

# intersect with bed file
subject <- with(tmb_bed_file, GRanges(chr, IRanges(start = start, end = end)))
query <- with(pbta_maf, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
pbta_tmb <- findOverlaps(query = query, subject = subject, type = "within")
pbta_tmb <- data.frame(pbta_maf[queryHits(pbta_tmb),], tmb_bed_file[subjectHits(pbta_tmb),])

# mutations per sample
pbta_tmb <- pbta_tmb %>%
  mutate(sample_name =  Tumor_Sample_Barcode) %>%
  group_by(sample_name) %>%
  mutate(num_var = n()) %>%
  mutate(tmb = num_var/77.46) %>%
  dplyr::select(sample_name, tmb) %>%
  unique()

# add short histology
pbta_histology <- read.delim(file.path(ref_dir, 'pbta', 'pbta-histologies.tsv'))
pbta_histology <- pbta_histology %>%
  dplyr::select(Kids_First_Biospecimen_ID, short_histology) %>%
  unique()
pbta_tmb <- pbta_tmb %>%
  inner_join(pbta_histology, by = c('sample_name'= 'Kids_First_Biospecimen_ID')) %>%
  dplyr::select(short_histology, sample_name, tmb)
colnames(pbta_tmb) <- c("Diseasetype", "Samplename", "TMBscore")

# update file with new filters
write.table(pbta_tmb, file = file.path(ref_dir, 'pbta-TMBscores_withdiseastype.txt'), quote = F, sep = "\t", row.names = F)

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
