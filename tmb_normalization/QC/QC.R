##### QC for the current run 
current_result<- readr::read_tsv("../results/TCGA_not_in_pbta_diseasetypes_and_samples_TMBscores.txt")
previous_result> previous_result <- readr::read_tsv("../references/TCGA_diseasetypes_and_samples_TMBscores.txt")
previous_result <- previous_result %>% rename(TMBscore_previous = TMBscore) %>% 
  mutate(BedLength_default = 77461748)

tcga_annotate <- readr::read_tsv("../results/tcga_not_in_pbta_bam_manifest.tsv")
tcga_bed_select <- readr::read_tsv("../results/tcga_not_in_pbta_bed_selected.tsv")
tcga_annotated <- tcga_annotate %>% left_join(tcga_bed_select) %>% rename(Samplename = associated_entities.entity_submitter_id) %>%
  select(Samplename, bed_selected)

combined <- left_join(current_result, previous_result)  %>% 
  mutate(expected_TMBscore = TMBscore_previous*BedLength_default / BedLength) %>%
  mutate(ratio = TMBscore / expected_TMBscore) %>% 
  left_join(tcga_annotated)

ggplot(combined, aes(bed_selected, ratio), fill = Diseasetype) + geom_boxplot() + geom_hline(aes(yintercept = 1), color = "red") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(combined, aes(Diseasetype, ratio), fill = Diseasetype) + geom_boxplot() + geom_hline(aes(yintercept = 1), color = "red") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

