####################################
# Tabulate RNA-Seq Pathway Analysis
####################################

# Function to return all results from RNA-Seq Analysis
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(optparse))

# parse params
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Directory e.g. data/PNOC008-04"),
  make_option(c("-o", "--output"), type = "character",
              help = "output excel filename i.e. PNOC008-04_summary"),
  make_option(c("-t", "--type"), type = "character",
              help = "text or excel")
)
opt <- parse_args(OptionParser(option_list = option_list))
topDir <- opt$input
fname <- opt$output
type <- opt$type
pnoc008_patient <- gsub('.*/','',topDir)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
gsea.dir <- file.path(ref_dir, 'gsea')

# output filename
fname <- file.path(topDir, "output", fname)

# read output from gsea enrichment
# gtex brain
pnoc008_vs_gtex_brain <- readRDS(file.path(gsea.dir, 'pnoc008_vs_gtex_brain.rds'))
pnoc008_vs_gtex_brain <- pnoc008_vs_gtex_brain[[pnoc008_patient]]

# pbta hgg
pnoc008_vs_pbta_hgg <- readRDS(file.path(gsea.dir, 'pnoc008_vs_pbta_hgg.rds'))
pnoc008_vs_pbta_hgg <- pnoc008_vs_pbta_hgg[[pnoc008_patient]]

# pbta all histologies
pnoc008_vs_pbta <- readRDS(file.path(gsea.dir, 'pnoc008_vs_pbta.rds'))
pnoc008_vs_pbta <- pnoc008_vs_pbta[[pnoc008_patient]]

# up/down pathways
pathway_df <- rbind(pnoc008_vs_gtex_brain$pathways, pnoc008_vs_pbta_hgg$pathways, pnoc008_vs_pbta$pathways)
pathway_df <- pathway_df %>%
  group_by(pathway, direction) %>%
  mutate(Freq = n()) %>%
  as.data.frame()
pathway_df.up <- pathway_df %>%
  filter(direction == "up") %>%
  as.data.frame()
pathway_df.down <- pathway_df %>%
  filter(direction == "down") %>%
  as.data.frame()

# up/down genes
genes_df <- rbind(pnoc008_vs_gtex_brain$genes, pnoc008_vs_pbta_hgg$genes, pnoc008_vs_pbta$genes)
genes_df <- genes_df %>%
  group_by(genes, diff_expr) %>%
  mutate(Freq = n()) %>%
  as.data.frame()
genes_df.up <- genes_df %>%
  filter(diff_expr == "up") %>%
  as.data.frame()
genes_df.down <- genes_df %>%
  filter(diff_expr == "down") %>%
  as.data.frame()

# write out to excel
# instead of excel, write output to four text files
if(type == 'excel'){
  fname <- paste0(fname, '.xlsx')
  write.xlsx(x = pathway_df.up, file = fname, sheetName = "Pathways_Up", row.names = F)
  write.xlsx(x = genes_df.up, file = fname, sheetName = "DE_Genes_Up", row.names = F, append = TRUE)
  write.xlsx(x = pathway_df.down, file = fname, sheetName = "Pathways_Down", row.names = F, append = TRUE)
  write.xlsx(x = genes_df.down, file = fname, sheetName = "DE_Genes_Down", row.names = F, append = TRUE)
} else {
  write.table(x = pathway_df.up, file = paste0(fname, "_Pathways_Up.txt"), row.names = F, quote = F, sep = "\t")
  write.table(x = genes_df.up, file = paste0(fname, "_DE_Genes_Up.txt"), row.names = F, quote = F, sep = "\t")
  write.table(x = pathway_df.down, file = paste0(fname, "_Pathways_Down.txt"), row.names = F, quote = F, sep = "\t")
  write.table(x = genes_df.down, file = paste0(fname, "_DE_Genes_Down.txt"), row.names = F, quote = F, sep = "\t")
}
