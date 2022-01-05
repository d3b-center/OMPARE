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
  print(module_of_interest)
  
  #### SubNetwork Generation ---------------------------------------------------------------
  
  # subset to module of interest
  subnetwork_each <- subnetwork %>% 
    dplyr::filter(Module == module_of_interest) 
  # Induce subnetwork
  subnetwork_graphed <-subnetwork_each %>%
    dplyr::select(Gene1, Gene2) %>% as.matrix() %>%
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
  list_of_qSigs <- Filter(nrow, list_of_qSigs) # remove subnet with no rows
  
  # Generate an ordered list of unique drugs (ordered by their WTCS scores)
  lapply(list_of_qSigs, function(x){
    drug_list <- x %>% 
      arrange(WTCS) %>%
      pull(Drug_Name) %>% unique() 
    nDrug <- length(drug_list)
    # Give all the drugs weighted score from 1-2 based on their ranking 
    weight_score <- (1 + (1.0 - (c(1:nDrug)/nDrug)))
    drug_score_list <- data.frame(drug_list, weight_score)
    
    # calculate closeness and betweenness for the module's subnetwork
    vSet0 <- V(subnetwork_graphed)$name
    sc1 <- closeness(subnetwork_graphed, vSet0)  
    sc2 <- betweenness(subnetwork_graphed, vSet0, directed=F) 
    sc3 <- page.rank(subnetwork_graphed, algo="prpack", vids=vSet0)$vector ### for normalization
    
    # define an empty list to store the info 
    scores <- list()
    # get closeness and betweenness for all drugs in the subnetwork
    for(m in 1:nrow(drug_score_list)){
      # iterate through all drugs
      drug_of_interest <- drug_score_list[m,1]
      
      # get the targets for this particular drug
      drug_targets <- x %>% 
        dplyr::filter(Drug_Name == drug_of_interest) %>%
        pull(hgnc_symbol) %>% unique() 
      
      # calculate closeness and betweenness for all drugs
      sc11 <- closeness(subnetwork_graphed, drug_targets)  
      sc12 <- betweenness(subnetwork_graphed, drug_targets, directed=F) 
      sc13 <- page.rank(subnetwork_graphed, algo="prpack", vids=drug_targets)$vector
      
      # generate a list 
      scores_each <- list(sc11, sc12, sc13)
      names(scores_each) <- c(paste0(drug_of_interest, "_sc11"), 
                              paste0(drug_of_interest, "_sc12"),
                              paste0(drug_of_interest, "_sc13"))
      scores <- c(scores, scores_each)
    }
    
    # Initiate data frame that would be used to store scores
    sScore_df <- data.frame(matrix(module_of_interest, nDrug*(nDrug-1)/2, 4))
    
    k <- 0
    for (i in 1:(nDrug-1)){
      # drug of interest
      drug_of_interest1 <- drug_list[i]
      # extract the weight of scores as well
      weight_of_drug1 <- drug_score_list %>% 
        dplyr::filter(drug_list == drug_of_interest1) %>%
        pull(weight_score) %>% as.numeric()
      
      for(j in (i+1):nDrug){
        k <- k+1
        # drug of interest
        drug_of_interest2 <- drug_list[j]
        
        # extract the weight of scores as well
        weight_of_drug2 <- drug_score_list %>% 
          dplyr::filter(drug_list == drug_of_interest2) %>%
          pull(weight_score) %>% as.numeric()
        
        # output the drug name to the synergy score table
        sScore_df[k,1] <- drug_of_interest1
        sScore_df[k,2] <- drug_of_interest2
        
        # extract closeness, betweeness and rank for drugs of interest 
        sc11 <- scores[[paste0(drug_of_interest1, "_sc11")]]
        sc12 <- scores[[paste0(drug_of_interest1, "_sc12")]]
        sc13 <- scores[[paste0(drug_of_interest1, "_sc13")]]
        
        # extract closeness, betweeness and rank for drugs of interest 
        sc21 <- scores[[paste0(drug_of_interest2, "_sc11")]]
        sc22 <- scores[[paste0(drug_of_interest2, "_sc12")]]
        sc23 <- scores[[paste0(drug_of_interest2, "_sc13")]]
        
        # calculate synergy score
        sc11 <- (sc11-min(sc1))/(max(sc1)-min(sc1)) 
        sc12 <- (sc12-min(sc2))/(max(sc2)-min(sc2)) 
        sc13 <- (sc13-min(sc3))/(max(sc3)-min(sc3))
        sc21 <- (sc21-min(sc1))/(max(sc1)-min(sc1)) 
        sc22 <- (sc22-min(sc2))/(max(sc2)-min(sc2)) 
        sc23 <- (sc23-min(sc3))/(max(sc3)-min(sc3))
        sctar1 <- (sum(sc11)+sum(sc12)+sum(sc13))/3.0
        sctar2 <- (sum(sc21)+sum(sc22)+sum(sc23))/3.0
        sScore <- sctar1 + sctar2
        
        # calculate the synergy scores
        vt0 <- sScore * weight_of_drug1 * weight_of_drug2
        sScore_df[k,3] <- vt0
      }
    }
    
    # write out the results 
    comparison_name <- x %>% 
      dplyr::select(comparison) %>% pull() %>% unique()
    colnames(sScore_df) <-c("drug1", "drug2", "synergy_score", "module")
    sScore_df <- sScore_df %>% arrange(desc(synergy_score)) %>%
      mutate(comparison = comparison_name)
    return(sScore_df)
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



