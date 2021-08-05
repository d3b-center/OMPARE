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

ggplot(combined, aes(bed_selected, ratio)) + geom_boxplot() + geom_hline(aes(yintercept = 1), color = "red") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(combined, aes(Diseasetype, ratio)) + geom_boxplot() + geom_hline(aes(yintercept = 1), color = "red") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##### A big difference between adult and kids were seen, therefore, compare the pediatric results to previous one to see 
# whether that is expected 

current_pediatric_results <- readr::read_tsv("../results/PBTA-TMBscores_withdiseastype.txt")
in_use_pediatric_results <- readr::read_tsv("../references/pbta-TMBscores_withdiseastype.txt") %>% 
  dplyr::rename(TMBscore_in_use = TMBscore)
old_corrected_pediatric_results <- readr::read_tsv("../references/pbta-corrected-TMBscores.txt") %>%
  dplyr::rename(TMBscore_old_corrected = TMB) %>% 
  dplyr::select(Samplename, experimental_strategy, TMBscore_old_corrected)

combined<- dplyr::left_join(current_pediatric_results, in_use_pediatric_results) %>% 
  dplyr::left_join(old_corrected_pediatric_results)

histologies <- readr::read_tsv("../references/pbta-histologies.tsv") %>% 
  dplyr::select(Kids_First_Biospecimen_ID, cohort, experimental_strategy) %>% 
  dplyr::rename(Samplename = Kids_First_Biospecimen_ID)

# calculate the ratio for different calling methods
combined <- combined %>% mutate(corrected_ratio = TMBscore / TMBscore_old_corrected) %>%
  mutate(in_use_vs_now_ratio = TMBscore_in_use / TMBscore) %>% dplyr::left_join(histologies) %>%
  # adding back the bed selected for normalization 
  dplyr::mutate(bed_selected = case_when(
    experimental_strategy == "WGS" ~ "CCDS.bed",
    experimental_strategy == "WXS" ~ "Strexome_targets_intersect_sorted_padded100.GRCh38.withCCDS.bed",
    experimental_strategy == "Targeted Sequencing" & cohort == "PNOC" ~ "Strexome_targets_intersect_sorted_padded100.GRCh38.withCCDS.bed",
    TRUE ~ "xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.withCCDS.bed"
  )) %>% filter(ExpStrategy != "Targeted Sequencing")


ggplot(combined, aes(Diseasetype, corrected_ratio)) + 
  geom_boxplot() + geom_hline(aes(yintercept = 1), color = "red") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(combined, aes(Diseasetype, in_use_vs_now_ratio)) + 
  geom_boxplot() + geom_hline(aes(yintercept = 1), color = "red") + 
  geom_boxplot() + geom_hline(aes(yintercept = 10), color = "blue") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(combined, aes(bed_selected, in_use_vs_now_ratio)) + 
  geom_boxplot() + geom_hline(aes(yintercept = 1), color = "red") + 
  geom_boxplot() + geom_hline(aes(yintercept = 10), color = "blue") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# calculate the expected TMB just due to the bed file size and see whether that is expected 
combined <- combined %>% mutate(expected_TMB = TMBscore_in_use * 77461748 / BedLength) %>% 
  mutate(loss_ratio = TMBscore / expected_TMB)

ggplot(combined, aes(bed_selected, loss_ratio)) + 
  geom_boxplot() + geom_hline(aes(yintercept = 1), color = "red")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(combined, aes(Diseasetype, loss_ratio)) + 
  geom_boxplot() + geom_hline(aes(yintercept = 1), color = "red")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

readr::write_tsv(combined, "../QC/combined_PBTA_TMBscore_comparison.tsv")
