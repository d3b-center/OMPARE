# Author: Run Jin
#
# Obtain Drug Target of Interest for each qSig output (comparing to GTEx, PBTA all or PBTA HGG)

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("igraph"))

#### Define Directories ----------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

analysis_dir <- file.path(root_dir, "code", "drug_synergy")
results_dir <- file.path(analysis_dir, "results", "synergy_score")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}

source(file.path(analysis_dir, "utils/get_synergy_score.R"))

#### Parse command line options ------------------------------------------------

option_list <- list(
  make_option(c("-g","--gtex_mapped"),type="character",
              help="Output file for drug-mapped gtex qSig subsetted subnetworks (.tsv) "),
  make_option(c("-p","--pbta_mapped"),type="character",
              help="Output file for drug-mapped pbta qSig subsetted subnetworks (.tsv) "),
  make_option(c("-h","--pbta_hgg_mapped"),type="character",
              help="Output file for drug-mapped pbta hgg subsetted subnetworks (.tsv) "),
  make_option(c("-f","--subnetwork"),type="character",
               help="File for subnetworks of module of interest (.tsv) "),
  make_option(c("-b","--output_path"),type="character",
              help="Path for results output ")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
output_path <- opt$output_path 

#### Read in files necessary for analyses -----------------------------------
gtex_qSig_subnet_mapped <- readr::read_tsv(opt$gtex_mapped) %>% mutate(comparison = "gtex_qSig")
pbta_qSig_subnet_mapped <- readr::read_tsv(opt$pbta_mapped) %>% mutate(comparison = "pbta_qSig")
pbta_hgg_qSig_subnet_mapped <- readr::read_tsv(opt$pbta_hgg_mapped) %>% mutate(comparison = "pbta_hgg_qSig")

subnetwork <- readr::read_tsv(opt$subnetwork)

#### Synergy Score calculated ----------------------------------------------------------
# following steps adapted from SynGeNet Publication [@doi:10.1038/s41540-019-0085-4]
# Supplementary information downloadable as `SynGeNet.R`

#### SubNetwork Generation ---------------------------------------------------------------
subnetwork_graphed <- subnetwork %>% select(Gene1, Gene2) %>% as.matrix() 

  # Induce background network
subnetwork_graphed <- graph.edgelist(subnetwork_graphed)

# Find the nodes of subnetwork (gene symbols in the subnetworks)
gene1_names <- subnetwork %>% pull(Gene1) %>% unique()
gene2_names <- subnetwork %>% pull(Gene2) %>% unique()
all_nodes <- union(gene1_names, gene2_names)

#### Synergy Score Calculation ---------------------------------------------------------------
list_of_qSigs <- list(gtex_qSig_subnet_mapped, pbta_qSig_subnet_mapped, pbta_hgg_qSig_subnet_mapped)

# Generate an ordered list of unique drugs (ordered by their WTCS scores)
lapply(list_of_qSigs, function(x){
  drug_list <- x %>% 
    arrange(WTCS) %>%
    pull(Drug_Name) %>% unique() 
  nDrug <- length(drug_list)
  
  # Give all the drugs weighted score from 1-2 based on their ranking 
  weight_score <- (1 + (1.0 - (c(1:nDrug)/nDrug)))
  drug_score_list <- data.frame(drug_list, weight_score)
  
  # Initiate dataframe that would be used to store scores
  sScore <- data.frame(matrix('test', nDrug*(nDrug-1)/2, 3))
  
  k <- 0
  targets <- x$hgnc_symbol
  drugs <- x$Drug_Name
  
  for (i in 1:(nDrug-1)){
    print(i)
    
    # drug of interest
    drug_of_interest1 <- drug_list[i]
    
    # gene targets of the drug
    drug_targets1 <- x %>% 
      dplyr::filter(Drug_Name == drug_of_interest1) %>%
      pull(hgnc_symbol) %>% unique() 
    
    # only keep targets that are in the subnetwork
    targets_in_sub1 <- intersect(drug_targets1, all_nodes)
    
    # subtract the weight of scores as well
    weight_of_drug1 <- drug_score_list %>% 
      dplyr::filter(drug_list == drug_of_interest1) %>%
      pull(weight_score) %>% as.numeric()
    
    for (j in (i+1):nDrug){
      k <- k+1
      # drug of interest
      drug_of_interest2 <- drug_list[j]
      
      # gene targets of the drug
      drug_targets2 <- x %>% 
        dplyr::filter(Drug_Name == drug_of_interest2) %>%
        pull(hgnc_symbol) %>% unique() 
      
      # only keep targets that are in the subnetwork
      targets_in_sub2 <- intersect(drug_targets2, all_nodes)
      
      # subtract the weight of scores as well
      weight_of_drug2 <- drug_score_list %>% 
        dplyr::filter(drug_list == drug_of_interest2) %>%
        pull(weight_score) %>% as.numeric()
      
      sScore[k,1] <- drug_of_interest1
      sScore[k,2] <- drug_of_interest2
      
      vt0 <- getSynScore2(targets_in_sub1, targets_in_sub2, subnetwork_graphed) * weight_of_drug1 * weight_of_drug2
      sScore[k,3] <- vt0
    
    }
  }
  
  # write out the results 
  comparison_name <- x %>% select(comparison) %>% pull() %>% unique()
  colnames(sScore) <-c("drug1", "drug2", "synergy_score")
  sScore %>% arrange(desc(synergy_score)) %>%
    readr::write_tsv(file.path(output_path,paste0(comparison_name, "_synergy_score.tsv")))
  })
    


