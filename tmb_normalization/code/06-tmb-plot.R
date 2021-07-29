# this script is adapted from https://github.com/d3b-center/OMPARE/
# blob/master/code/patient_level_analyses/utils/tmb_profile.R in the OMPARE repo

library(tidyverse)
library(GenomicRanges)

# Load the required files 
ped_tmb <- data.table::fread('../results/PBTA-TMBscores_withdiseastype.txt') %>% 
  select(Diseasetype, Samplename, TMBscore)
adult_tmb <- data.table::fread('../results/TCGA_not_in_pbta_diseasetypes_and_samples_TMBscores.txt') %>%
  select(Diseasetype, Samplename, TMBscore)
adult_tmb_in_pbta <- data.table::fread('../results/TCGA_in_pbta_diseasetypes_and_samples_TMBscores.txt') %>%
  select(Diseasetype, Samplename, TMBscore)

tmb_bed_file <- data.table::fread("../references/ashion_confidential_exome_v2_2nt_pad.Gh38.bed")
tmb_bed_file <- tmb_bed_file[,1:3]
colnames(tmb_bed_file)  <- c('chr', 'start', 'end')
tmb_bed_file$chr <- sub("^", "chr", tmb_bed_file$chr )

# calculate the length of the bed file
bed_length<-0
for (i in 1:nrow(tmb_bed_file)) {
  distance <- as.numeric(tmb_bed_file[i,3]) - as.numeric(tmb_bed_file[i,2])
  bed_length <- bed_length + distance
}
  
## read in parameters
var_class = c('Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins')
vaf_cutoff = 0.05
var_count = 3
tumor_depth = 25

# read mutect2 for TMB profile

myMutData_32 <- data.table::fread("../references/pnoc008_32.mutect2_somatic.norm.annot.protected.maf",skip = 1, stringsAsFactors = F)
myMutData_33 <- data.table::fread("../references/pnoc008_33.mutect2_somatic.norm.annot.protected.maf",skip = 1, stringsAsFactors = F)
myMutData_34 <- data.table::fread("../references/pnoc008_34.mutect2_somatic.norm.annot.protected.maf",skip = 1, stringsAsFactors = F)

mutFiles <- list(myMutData_32, myMutData_33, myMutData_34)

# apply filters using Friends of Cancer Research Project Standards 
myMutData_list <- lapply(mutFiles, function(x){
  x %>%
  dplyr::mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
    dplyr::filter(Variant_Classification %in% var_class,
         t_depth >= tumor_depth,
         vaf >= vaf_cutoff,
         t_alt_count >= var_count) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
})

# intersect with bed file
TMB_list <- lapply(myMutData_list, function(x){
  myMutData=x
  subject <- with(tmb_bed_file, GRanges(chr, IRanges(start = start, end = end)))
  query <- with(myMutData, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
  res <- findOverlaps(query = query, subject = subject, type = "within")
  res <- data.frame(myMutData[queryHits(res),], tmb_bed_file[subjectHits(res),])
  TMB_each <- (nrow(res))*1000000/bed_length
  return(TMB_each)
  }
)
  
  ped_tmb$Type <- "Pediatric"
  adult_tmb$Type <- "Adult"
  adult_tmb_in_pbta$Type <- "Adult in PBTA"
  
  tmbScores <- rbind(ped_tmb, adult_tmb, adult_tmb_in_pbta)
  
  # count median per histology
  disease.order <- tmbScores %>%
    group_by(Type, Diseasetype) %>%
    summarise(median = median(TMBscore), count = n()) %>%
    filter(count > 2) %>%
    arrange(desc(median)) %>%
    .$Diseasetype
  
  # subset to disease type with > 2 samples
  tmbScores <- tmbScores %>%
    filter(Diseasetype %in% disease.order)
  tmbScores$Diseasetype <- factor(tmbScores$Diseasetype, levels = disease.order)
  
  # Plot it
  ggplot_list <- lapply(TMB_list, function(x){
  ggplot(tmbScores, aes(Diseasetype, TMBscore, fill = Type)) + 
    geom_boxplot() + theme_bw() + scale_y_log10(breaks = c(.25, 1, 10, 100, 500)) + 
    scale_fill_manual(values = c("blue", "orange", "green")) + 
    xlab("Disease") + ylab("Mutations per MB") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    geom_hline(yintercept = x, linetype = 2, color = 'gray30') +
    annotate("text", x = 35, y = max(tmbScores$TMBscore) - 50, 
             label = "- - - Patient TMB", size = 3, 
             fontface = 'italic', color = "gray30")
  })

  ggsave("../plots/PNOC008-32-TMB-plot.png", ggplot_list[[1]])
  ggsave("../plots/PNOC008-33-TMB-plot.png", ggplot_list[[2]])
  ggsave("../plots/PNOC008-34-TMB-plot.png", ggplot_list[[3]])
