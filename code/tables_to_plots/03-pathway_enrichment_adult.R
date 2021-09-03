suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient identifier, i.e. PNOC008-XX"),
  make_option(c("--output_dir"), type = "character",
              help = "Output directory")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
pnoc008_sample_of_interest <- opt$patient
output_dir <- opt$output_dir

## P6 pathway enrichment
### Read in the files and arrange for plotting
pathway_analysis_adult <- readRDS(file.path(output_dir, "pathway_analysis_adult.rds"))

shared_pathway_adult <- pathway_analysis_adult$shared_pathways %>%
  # arrange by counts
  arrange(desc(Sample.count.per.pathway))  %>% 
  select(pathway, Sample.count.per.pathway, direction) %>% 
  distinct() %>%
  # take the top 10 for each group
  group_by(direction) %>% 
  dplyr::slice(1:10) %>% 
  arrange(Sample.count.per.pathway) %>% 
  mutate(direction = factor(direction, levels = c("up", "down")),
         pathway = factor(pathway, levels = unique(pathway)))

### Plot for pathway enrichment
p <- ggplot(shared_pathway_adult, aes(pathway, y = Sample.count.per.pathway, fill = direction)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() +
    xlab("") + 
    ylab("Count of Enriched Pathways in 20 Transcriptomically Similar Patients") + 
    theme(plot.margin = unit(c(1, 5, 1, 7), "cm")) + 
    scale_fill_manual(name = "Direction", 
                      values = c("up" = "red", "down" = "forest green")) + 
    ggtitle(paste0("Comparison against adult")) 

# save in patient of interest's output folder
ggsave(plot = p, 
       filename = file.path(output_dir, "pathway_analysis_adult.pdf"), 
       height = 6, width = 15)
