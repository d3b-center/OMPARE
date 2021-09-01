suppressPackageStartupMessages({
  library(viridis)
  library(tidyverse)
  library(ggplot2)
  library(optparse)
})

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source theme
source(file.path(patient_level_analyses_utils, "pubTheme.R"))

# Parse command line options
option_list <- list(
  make_option(c("--combined_synergy"),type="character",
              help="Path and file name for all combined synergy score (.tsv)"),
  make_option(c("--output_file"),type="character",
              help="Path and file name for output bubble plot (.pdf)")
)
opt <- parse_args(OptionParser(option_list = option_list))
combined_synergy <- opt$combined_synergy
output_file <- opt$output_file

# read output of  
all_combined <- read.delim('results/PNOC008-38/output/drug_synergy/combined_qSig_synergy_score.tsv')
all_combined <- read.delim(combined_synergy)

# create drug pairs
all_combined <- all_combined %>%
  mutate(drug_pair = paste(drug1, drug2, sep = "_")) %>%
  dplyr::select(drug_pair, module, synergy_score, comparison) %>%
  unique()
  
# modify comparison names
all_combined$comparison[all_combined$comparison == "pbta_qSig"] <- "PBTA"
all_combined$comparison[all_combined$comparison == "pbta_hgg_qSig"] <- "PBTA (HGG)"
all_combined$comparison[all_combined$comparison == "gtex_qSig"] <- "GTEx"

# pick top 10 drug pairs per comparator 
test <- all_combined %>%
  group_by(comparison) %>%
  arrange(desc(synergy_score)) %>%
  slice(1:10)

# sort by alphabetical order
test$drug_pair <- factor(test$drug_pair, levels = unique(sort(test$drug_pair)))

# plot
p <- ggplot(data = test, aes(x = comparison, y = drug_pair, fill = comparison)) +
  geom_point(aes(size = synergy_score),
             shape = 21, 
             color = "black", alpha = 0.8) +
  scale_fill_manual(values = viridis_pal(option = "plasma")(4)) + 
  scale_size(range = c(.1, 15)) + theme_Publication() +
  guides(fill = "none") + xlab("") + ylab("") + 
  ggtitle("Drug pair synergy scores")
ggsave(filename = output_file, width = 10, height = 8)

