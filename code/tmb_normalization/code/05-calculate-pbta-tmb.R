library(tidyverse)
library(GenomicRanges)

# filters
var_class = c('Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',  'In_Frame_Del', 'In_Frame_Ins')
vaf_cutoff = 0.05
var_count = 3
tumor_depth = 25

# set root directory and other directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
results_dir <- file.path(root_dir, 'tmb_normalization', 'results')
ref_dir <- file.path(root_dir, 'tmb_normalization', 'references')

# read in the PBTA mutect2 file 
pbta_mutect2 <- data.table::fread(file.path(ref_dir, 'pbta-snv-mutect2.vep.maf.gz'))

# read in the histologies file
histologies <- read.delim(file.path(ref_dir, "pbta-histologies.tsv"), header = T, sep = "\t", stringsAsFactors = F)

# Take a look at all the samples that are in there 
pbta_mutect2_samples <- pbta_mutect2$Tumor_Sample_Barcode %>% unique()

# Filter the histologies file to contain only the specimens that have mutect2 results
pbta_mutect2_samples_histology <- histologies %>% 
  dplyr::filter(Kids_First_Biospecimen_ID %in% pbta_mutect2_samples)

# See what are the cohorts and experimental strategy
pbta_mutect2_samples_histology$cohort %>% unique()
pbta_mutect2_samples_histology$experimental_strategy %>% unique()

# read in the default bed file for all PBTA samples
pbta_bed <- read.delim(file.path(ref_dir, "pbta_bed_files", "xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed"))
pbta_bed <- pbta_bed[,1:3]
colnames(pbta_bed)  <- c('chr', 'start', 'end')

# change the histology file to contain common column as pbta_maf
pbta_mutect2_samples_histology <- pbta_mutect2_samples_histology %>%
  dplyr::mutate(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Tumor_Sample_Barcode, short_histology, experimental_strategy, cohort)

# filter pbta maf file 
pbta_maf <- pbta_mutect2 %>%
  group_by(Tumor_Sample_Barcode) %>%
  mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
  filter(Variant_Classification %in% var_class,
         t_depth >= tumor_depth,
         vaf >= vaf_cutoff,
         t_alt_count >= var_count) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)

# annotate maf with relevant histologies information 
pbta_maf <- pbta_maf %>% left_join(pbta_mutect2_samples_histology)

# calculate bed length for the default bed file 
bed_length<-0
for (i in 1:nrow(pbta_bed)) {
  distance <- as.numeric(pbta_bed[i,3]) - as.numeric(pbta_bed[i,2])
  bed_length = bed_length + distance
}

# calculate maf result using the default bed file 
# intersect with bed file
subject <- with(pbta_bed, GRanges(chr, IRanges(start = start, end = end)))
query <- with(pbta_maf, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
pbta_tmb <- findOverlaps(query = query, subject = subject, type = "within")
pbta_tmb <- data.frame(pbta_maf[queryHits(pbta_tmb),], pbta_bed[subjectHits(pbta_tmb),])

# mutations per sample
pbta_tmb <- pbta_tmb %>%
  mutate(sample_name =  Tumor_Sample_Barcode) %>%
  group_by(sample_name) %>%
  mutate(num_var = n()) %>%
  mutate(tmb = num_var*1000000/bed_length) %>%
  dplyr::select(sample_name, short_histology, cohort, experimental_strategy, tmb) %>%
  unique() %>% 
  dplyr::mutate(BedLength = bed_length)

colnames(pbta_tmb) <- c("Samplename", "Diseasetype", "Cohort", "ExpStrategy", "TMBscore", "BedLength")

# update file with new filters and save to data directory
write.table(pbta_tmb, file = file.path(data_dir, 'tmb', 'PBTA-TMBscores_withdiseastype.txt'), quote = F, sep = "\t", row.names = F)

