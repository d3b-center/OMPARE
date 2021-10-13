# Author: Run Jin
# Obtain Drug Target of Interest for each qSig output (comparing to GTEx, PBTA all or PBTA HGG)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(igraph)
})

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# module path
module_dir <- file.path(root_dir, "code", "drug_synergy")
source(file.path(module_dir, "utils", "get_synergy_score.R"))

# Parse command line options
option_list <- list(
  make_option(c("-g","--gtex_mapped"),type="character",
              help="Output file for drug-mapped gtex qSig subsetted subnetworks (.tsv) "),
  make_option(c("-p","--pbta_mapped"),type="character",
              help="Output file for drug-mapped pbta qSig subsetted subnetworks (.tsv) "),
  make_option(c("-h","--pbta_hgg_mapped"),type="character",
              help="Output file for drug-mapped pbta hgg subsetted subnetworks (.tsv) "),
  make_option(c("-f","--subnetwork"),type="character",
               help="File for subnetworks of module of interest (.tsv) "),
  make_option(c("-m","--subnetwork_mapped"),type="character",
              help="File for subnetworks with mapped drug information of module of interest (.tsv) "),
  make_option(c("-b","--output_gtex"),type="character",
              help="Path and file name for gtex synergy score (.tsv) "), 
  make_option(c("-c","--output_pbta"),type="character",
              help="Path and file name for pbta synergy score (.tsv) "), 
  make_option(c("-d","--output_pbta_hgg"),type="character",
              help="Path and file name for pbta hgg synergy score (.tsv) "),
  make_option(c("-e","--output_combined"),type="character",
              help="Path and file name for all combined synergy score (.tsv) ")
)
opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
output_path <- opt$output_path 

#### Read in files necessary for analyses -----------------------------------
gtex_qSig_subnet_mapped <- readr::read_tsv(opt$gtex_mapped) %>% mutate(comparison = "gtex_qSig")
pbta_qSig_subnet_mapped <- readr::read_tsv(opt$pbta_mapped) %>% mutate(comparison = "pbta_qSig")
pbta_hgg_qSig_subnet_mapped <- readr::read_tsv(opt$pbta_hgg_mapped) %>% mutate(comparison = "pbta_hgg_qSig")

if(nrow(gtex_qSig_subnet_mapped) == 0 |
   nrow(pbta_qSig_subnet_mapped) == 0 |
   nrow(pbta_hgg_qSig_subnet_mapped) == 0){
  print("No qSig output")
  stop()
}

subnetwork <- readr::read_tsv(opt$subnetwork)

#### Synergy Score calculated ----------------------------------------------------------
# following steps adapted from SynGeNet Publication [@doi:10.1038/s41540-019-0085-4]
# Supplementary information downloadable as `SynGeNet.R`

# Calculate SynergyScore for all subnetworks
module_list <- subnetwork %>% pull(Module) %>% unique()

with_module_list <- lapply(module_list, function(y){
  
  # iterate through each module
  module_of_interest <- y
  
  #### SubNetwork Generation ---------------------------------------------------------------
  
  # subset to module of interest
  subnetwork_each <- subnetwork %>% 
    filter(Module == module_of_interest) 
  # Induce subnetwork
  subnetwork_graphed <-subnetwork_each %>%
    select(Gene1, Gene2) %>% as.matrix() %>%
    graph.edgelist()
  
  # Find the nodes of subnetwork (gene symbols in the subnetworks)
  gene1_names <- subnetwork_each %>% pull(Gene1) %>% unique()
  gene2_names <- subnetwork_each %>% pull(Gene2) %>% unique()
  all_nodes <- union(gene1_names, gene2_names)
  
  # filter all three qSig files to module of interest as well
  gtex_qSig_subnet_mapped_each <- gtex_qSig_subnet_mapped %>% 
    filter(module == module_of_interest)
  pbta_qSig_subnet_mapped_each <- pbta_qSig_subnet_mapped %>% 
    filter(module == module_of_interest)
  pbta_hgg_qSig_subnet_mapped_each <- pbta_hgg_qSig_subnet_mapped %>% 
    filter(module == module_of_interest)
  
  #### Synergy Score Calculation ---------------------------------------------------------------
  list_of_qSigs <- list(gtex_qSig_subnet_mapped_each, pbta_qSig_subnet_mapped_each, pbta_hgg_qSig_subnet_mapped_each)
  
  # Generate an ordered list of unique drugs (ordered by their WTCS scores)
  lapply(list_of_qSigs, function(x){
    drug_list <- x %>% 
      arrange(WTCS) %>%
      pull(Drug_Name) %>% unique() 
    nDrug <- length(drug_list)
    # Give all the drugs weighted score from 1-2 based on their ranking 
    weight_score <- (1 + (1.0 - (c(1:nDrug)/nDrug)))
    drug_score_list <- data.frame(drug_list, weight_score)
    
    # Initiate data frame that would be used to store scores
    sScore <- data.frame(matrix(module_of_interest, nDrug*(nDrug-1)/2, 4))
    
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
      
      for(j in (i+1):nDrug){
        print(j)
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
    colnames(sScore) <-c("drug1", "drug2", "synergy_score", "module")
    sScore <- sScore %>% arrange(desc(synergy_score)) %>%
      mutate(comparison = comparison_name)
    return(sScore)
    })
  }
)  

each_module_combined <- lapply (with_module_list, function(x){
  do.call(rbind, x)
}) 

all_combined <- do.call(rbind, each_module_combined)

# annotate MOA to all the drugs
subnetwork_mapped <- readr::read_tsv(opt$subnetwork_mapped)

drug1_moa <- subnetwork_mapped %>% select(Drug_Name, MOA, hgnc_symbol) %>% distinct() %>% 
  rename(drug1 = Drug_Name, drug1_MOA=MOA, drug1_target_hgnc_symbol = hgnc_symbol)

drug2_moa <- subnetwork_mapped %>% select(Drug_Name, MOA, hgnc_symbol) %>% distinct() %>% 
  rename(drug2 = Drug_Name, drug2_MOA=MOA, drug2_target_hgnc_symbol = hgnc_symbol)


# writing out results 
all_combined %>% filter(comparison == "gtex_qSig") %>% 
  dplyr::left_join(drug1_moa) %>% 
  dplyr::left_join(drug2_moa) %>%
  arrange(desc(synergy_score)) %>% 
  readr::write_tsv(opt$output_gtex)

all_combined %>% filter(comparison == "pbta_qSig") %>%
  dplyr::left_join(drug1_moa) %>% 
  dplyr::left_join(drug2_moa) %>%
  arrange(desc(synergy_score)) %>%
  readr::write_tsv(opt$output_pbta)

all_combined %>% filter(comparison == "pbta_hgg_qSig")%>%
  dplyr::left_join(drug1_moa) %>% 
  dplyr::left_join(drug2_moa) %>%
  arrange(desc(synergy_score)) %>%
  readr::write_tsv(opt$output_pbta_hgg)

all_combined %>% arrange(desc(synergy_score)) %>%
  dplyr::left_join(drug1_moa) %>% 
  dplyr::left_join(drug2_moa) %>%
  readr::write_tsv(opt$output_combined)



