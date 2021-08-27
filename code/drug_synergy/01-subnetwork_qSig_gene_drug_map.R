# Author: Run Jin
#
# Obtain Gene Drug Mapping for Subnetworks of Patients of Interest

# BiocManager::install("drugTargetInteractions")
# BiocManager::install("biomaRt")

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("drugTargetInteractions"))
suppressPackageStartupMessages(library("biomaRt"))

#### Define Directories ----------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

ref_dir <- file.path(root_dir, "references")

results_dir <- file.path(root_dir, "code", "drug_synergy", "results", "drug_gene_map")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}

#### Parse command line options ------------------------------------------------

option_list <- list(
  make_option(c("-i", "--interaction"),type="character",
              help="Gene interaction file for all modules (.tsv) "),
  make_option(c("-r","--enrichment_nes"),type="character",
              help="cluster vs. modules enrichment score file (.tsv) "),
  make_option(c("-c","--cluster"),type="character",
              help="File with cluster assignments for samples (.tsv)"),
  make_option(c("-s","--sample_interest"),type="character",
              help="sample of interest to run this analysis "),
  make_option(c("-l","--chemblDb_path"),type="character",
              help="Path to chemblDb sqlite database"),
  make_option(c("-e", "--gtex_qSig"),type="character",
              help="qSig output from CEMiTools comparing to GTEx brain normal (.txt) "),
  make_option(c("-b","--pbta_qSig"),type="character",
              help="qSig output from CEMiTools comparing to all PBTA samples (.txt) "),
  make_option(c("-a","--pbta_hgg_qSig"),type="character",
              help="qSig output from CEMiTools comparing to PBTA HGG samples (.txt)"),
  make_option(c("-f","--subnetwork"),type="character",
              help="Output file for subnetworks of module of interest (.tsv) "),
  make_option(c("-o","--subnetwork_mapped"),type="character",
              help="Output file for drug-mapped subnetworks (.tsv) "),
  make_option(c("-g","--gtex_mapped"),type="character",
              help="Output file for drug-mapped gtex qSig subsetted subnetworks (.tsv) "),
  make_option(c("-p","--pbta_mapped"),type="character",
              help="Output file for drug-mapped pbta qSig subsetted subnetworks (.tsv) "),
  make_option(c("-d","--pbta_hgg_mapped"),type="character",
              help="Output file for drug-mapped pbta hgg subsetted subnetworks (.tsv) ")
)

opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
sample_interest <- opt$sample_interest

#### Read in files necessary for analyses -----------------------------------
interaction_df <- readr::read_tsv(opt$interaction)

enrichment_nes_df <- readr::read_tsv(opt$enrichment_nes) %>%
  dplyr::rename(cluster_1 = `1`,
                cluster_2 = `2`,
                cluster_3 = `3`)

cluster_df <- readr::read_tsv(opt$cluster)

gtex_qSig <- readr::read_tsv(opt$gtex_qSig)
pbta_qSig <- readr::read_tsv(opt$pbta_qSig)
pbta_hgg_qSig <- readr::read_tsv(opt$pbta_hgg_qSig)

#### Run drug annnotation for subnetwork -------------------------------------

# Find the cluster assignment of samples of interest 
cluster_assignment <- cluster_df %>% 
  dplyr::filter(Sample.Names == sample_interest ) %>%
  pull(CC)

cluster_assignment <- paste0("cluster_", cluster_assignment)
  
# Find the modules positive correlated with the cluster
enrichment_nes_df <- enrichment_nes_df %>%
  dplyr::select(pathway, all_of(cluster_assignment))
colnames(enrichment_nes_df) = c("modules", "cluster_of_interest")

module_selected <- enrichment_nes_df %>% 
  arrange(desc(cluster_of_interest)) %>%
  slice(1) %>%
  pull(modules)

# Use the module selected to subset the gene interaction file for all modules and generate sub-network
sub_network <- interaction_df %>%
  dplyr::filter(Module == module_selected)

write_tsv(sub_network, opt$subnetwork)

# Get the union of all genes in the subnetwork
gene1_list <- sub_network %>% pull(Gene1) %>% unique()
gene2_list <- sub_network %>% pull(Gene2) %>% unique()

all_genes_in_subnetwork <- union(gene1_list, gene2_list) %>% unique()

# use the vector of all genes in the subnetwork as input to getSymEnsUp function.
chemblDbPath <- opt$chemblDb_path
resultsPath <- system.file("extdata", "results", package = "drugTargetInteractions")
config <- genConfig(chemblDbPath = chemblDbPath, resultsPath = resultsPath)

# generate a mapping of uniprot ids and gene symbols only once save as rds
map_file <- file.path(ref_dir, "uniprot_genesymbol_map.rds")
if(!file.exists(map_file)){
  hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
  annotation.uni <- getBM(attributes=c("uniprotswissprot", "hgnc_symbol"), mart=hsmart)
  saveRDS(annotation.uni, file = map_file)
} else {
  annotation.uni <- readRDS(map_file)
}

# filter to genes of interest
annotation.uni <- annotation.uni %>%
  filter(hgnc_symbol %in% all_genes_in_subnetwork )
queryBy <- list(molType = "protein", idType = "UniProt_ID", ids = annotation.uni$uniprotswissprot)
qresult2 <- drugTargetAnnot(queryBy, config=config)
qresult2 <- qresult2 %>%
  filter(!is.na(First_Approval))

# add gene name back to table
qresult2 <- qresult2 %>%
  inner_join(annotation.uni, by = c("UniProt_ID" = "uniprotswissprot"))

write_tsv(qresult2, opt$subnetwork_mapped)

#### Intersect subnetwork and qSig output-------------------------------------

# First select drugs that are negatively correlated (WTCS<0) and keep the WTCS scores
gtex_qSig_df <- gtex_qSig %>% 
  dplyr::filter(WTCS<0) %>% mutate(Drug_Name = toupper(pert)) %>% 
  dplyr::select(Drug_Name, WTCS)
gtex_qSig_drugs <- gtex_qSig_df %>% pull(Drug_Name) %>% unique()

pbta_qSig_df <- pbta_qSig %>% 
  dplyr::filter(WTCS<0) %>% mutate(Drug_Name = toupper(pert)) %>% 
  dplyr::select(Drug_Name, WTCS)
pbta_qSig_drugs <- pbta_qSig_df %>% pull(Drug_Name) %>% unique()

pbta_hgg_qSig_df <- pbta_hgg_qSig %>% 
  dplyr::filter(WTCS<0) %>% mutate(Drug_Name = toupper(pert)) %>% 
  dplyr::select(Drug_Name, WTCS)
pbta_hgg_qSig_drugs <- pbta_hgg_qSig_df %>% pull(Drug_Name) %>% unique()

# Combine the results with WTCS scores and save the results
subnetwork_gtex_qSig <- qresult2 %>% 
  dplyr::filter(Drug_Name %in% gtex_qSig_drugs) %>%
  dplyr::left_join(gtex_qSig_df) %>% distinct() %>% 
  readr::write_tsv(opt$gtex_mapped)

subnetwork_pbta_qSig <- qresult2 %>% 
  dplyr::filter(Drug_Name %in% pbta_qSig_drugs) %>%
  dplyr::left_join(pbta_qSig_df) %>% distinct() %>% 
  readr::write_tsv(opt$pbta_mapped)

subnetwork_pbta_hgg_qSig <- qresult2 %>% 
  dplyr::filter(Drug_Name %in% pbta_hgg_qSig_drugs) %>%
  dplyr::left_join(pbta_hgg_qSig_df) %>% distinct() %>% 
  readr::write_tsv(opt$pbta_hgg_mapped)

