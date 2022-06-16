# function to create histology file
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

create_histologies_subset <- function(hist_file, group_filter, cohort_filter, output_dir, prefix){
  
  # create output directory
  dir.create(output_dir, showWarnings = F, recursive = T)
  
  # create short_histology for GTEx for comparisons
  # filter to specified cohort
  hist_file <- hist_file %>%
    dplyr::mutate(short_histology = ifelse(cohort == "GTEx", gtex_group, short_histology)) %>%
    filter(cohort %in% cohort_filter)
  
  # filter to specified group
  if(!is.null(group_filter)){
    hist_file <- hist_file %>%
      filter(short_histology %in% group_filter)
  }
  
  # for adult tumors, harmonize sample_id
  if(prefix == "adult_tumors"){
    hist_file <- hist_file %>%
      mutate(sample_id = gsub('[R]-[0-9A-Z]{4}-[0-9]{2}', '', Kids_First_Biospecimen_ID))
  }
    
  # save subset
  hist_output <- file.path(output_dir, paste(prefix, "histologies.tsv", sep = "_"))
  readr::write_tsv(hist_file, file = hist_output)
  
  return(hist_file)
}
