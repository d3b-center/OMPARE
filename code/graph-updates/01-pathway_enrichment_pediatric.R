suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient identifier, i.e. PNOC008-XX")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
pnoc008_sample_of_interest <- opt$patient

# Define directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
output_dir <- file.path(root_dir, "results", pnoc008_sample_of_interest, "output")

## P5 pathway enrichment
### Read in the files and arrange for plotting
pathway_analysis_pediatric <- readRDS(file.path(output_dir, "pathway_analysis_pediatric.rds"))
shared_pathway_pediatric <- pathway_analysis_pediatric$shared_pathways
shared_pathway_pediatric <- shared_pathway_pediatric %>% 
  filter(sample_name == pnoc008_sample_of_interest) %>%
  mutate(padj = as.numeric(padj),
         direction = factor(direction, levels = c("up", "down")),
         pathway = factor(pathway, levels = unique(pathway))) %>%
  arrange(padj) %>% 
  arrange(direction)

### Plot for pathway enrichment
p <- ggplot(shared_pathway_pediatric, aes(pathway, y = (-1)*log10(padj), fill = direction)) +
  geom_bar(stat="identity") + coord_flip() + theme_bw() +
  xlab("") + 
  ylab("-log10 Adj. P-Value") + 
  theme(plot.margin = unit(c(1, 5, 1, 7), "cm")) + 
  scale_fill_manual(name = "Direction", 
                    values = c("down" = "forest green", "up" = "red")) + 
  ggtitle(paste0("Comparison against pediatric")) 

# save in patient of interest's output folder
ggsave(plot = p, 
       filename = file.path(output_dir, "pathway_analysis_pediatric.pdf"),
       height = 6, width = 15)
