####################################
# Tabulate RNA-Seq Pathway Analysis
####################################

# Function to return all results from RNA-Seq Analysis
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(optparse))

# parse params
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Directory e.g. data/PNOC008-04"),
  make_option(c("-o", "--output"), type = "character",
              help = "output excel file with extension i.e. output.xlsx")
)
opt <- parse_args(OptionParser(option_list = option_list))
topDir <- opt$input
fname <- opt$output
pnoc008.patient <- gsub('.*/','',topDir)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
gsea.dir <- file.path(ref_dir, 'GSEA')

# output filename
fname <- file.path(topDir, "output", fname)

# read output from gsea enrichment
# gtex brain
GTExBrain <- readRDS(file.path(gsea.dir, 'PNOC008_vs_GTExBrain.RDS'))
GTExBrain <- GTExBrain[[pnoc008.patient]]

# pbta hgg
PBTA_HGG <- readRDS(file.path(gsea.dir, 'PNOC008_vs_PBTA_HGG.RDS'))
PBTA_HGG <- PBTA_HGG[[pnoc008.patient]]

# pbta all histologies
PBTA_All <- readRDS(file.path(gsea.dir, 'PNOC008_vs_PBTA.RDS'))
PBTA_All <- PBTA_All[[pnoc008.patient]]

# up/down pathways
pathway.df <- rbind(GTExBrain$pathways, PBTA_HGG$pathways, PBTA_All$pathways)
pathway.df <- pathway.df %>%
  group_by(Pathway, Direction) %>%
  mutate(Freq = n()) %>%
  as.data.frame()
pathway.df.up <- pathway.df %>%
  filter(Direction == "Up") %>%
  as.data.frame()
pathway.df.down <- pathway.df %>%
  filter(Direction == "Down") %>%
  as.data.frame()

# up/down genes
genes.df <- rbind(GTExBrain$genes, PBTA_HGG$genes, PBTA_All$genes)
genes.df <- genes.df %>%
  group_by(Gene_name, DE) %>%
  mutate(Freq = n()) %>%
  as.data.frame()
genes.df.up <- genes.df %>%
  filter(DE == "Up") %>%
  as.data.frame()
genes.df.down <- genes.df %>%
  filter(DE == "Down") %>%
  as.data.frame()

# write out to excel
write.xlsx(x = pathway.df.up, file = fname, sheetName = "Pathways_Up", row.names = F)
write.xlsx(x = genes.df.up, file = fname, sheetName = "DE_Genes_Up", row.names = F, append = TRUE)
write.xlsx(x = pathway.df.down, file = fname, sheetName = "Pathways_Down", row.names = F, append = TRUE)
write.xlsx(x = genes.df.down, file = fname, sheetName = "DE_Genes_Down", row.names = F, append = TRUE)
