# function: generate all findings table 
suppressPackageStartupMessages({
  library(tidyverse)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# read gencode from OpenPedCan
gencode_gtf <- rtracklayer::import(con = file.path(data_dir, 'OpenPedCan-analysis/data/gencode.v27.primary_assembly.annotation.gtf.gz'))
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  unique()

all_findings <- function(sample_info, annotated_maf, filtered_fusions, filtered_cnv, rnaseq_analysis_output, snv_caller) {
  
  # Add Mutations
  if(!is.null(annotated_maf)){
    # combine annotated maf with sample info
    annotated_maf <- annotated_maf %>%
      mutate(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
      inner_join(sample_info, by = 'Kids_First_Biospecimen_ID')
    
    # annotate mutations, add Aberration and Details
    tmpMut <- annotated_maf %>%
      dplyr::mutate(Aberration = Hugo_Symbol) %>%
      filter(HGVSp_Short != "") %>%
      dplyr::mutate(Details = paste0('Variant: ', Variant_Classification, " | HGVSp: ", HGVSp_Short)) 
    
    # only keep Pfam and SMART domains, collapse DOMAINS to string
    tmpMut <- tmpMut %>%
      dplyr::mutate(DOMAINS = ifelse(DOMAINS == "", NA, strsplit(as.character(DOMAINS), ","))) %>% 
      unnest(DOMAINS) %>%
      dplyr::mutate(match = str_detect(DOMAINS, 'Pfam_domain|SMART_domains')) %>%
      dplyr::mutate(DOMAINS = ifelse(match, DOMAINS, NA)) %>%
      group_by(Aberration, Type, Details) %>%
      dplyr::mutate(DOMAINS = toString(na.omit(DOMAINS))) %>%
      unique()
    
    # add COSMIC
    tmpMut <- tmpMut %>%
      dplyr::mutate(Existing_variation = str_detect(Existing_variation, 'COSM')) %>%
      dplyr::mutate(Existing_variation = ifelse(Existing_variation, "Cosmic_Variant", ""))
    
    # now add Variant Properties depending on snv_caller
    if(snv_caller == "all"){
      tmpMut <- tmpMut %>%
        dplyr::select(Aberration, Type, Details, SIFT, DOMAINS, Existing_variation, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id) %>%
        unique()
    } else {
      tmpMut <- tmpMut %>%
        dplyr::select(Aberration, Type, Details, t_alt_count, t_depth, SIFT, DOMAINS, Existing_variation, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id) %>%
        unique()
    }
    
    # collapse into Variant_Properties
    tmpMut <- reshape2::melt(as.data.frame(tmpMut), id.vars = c('Aberration', 'Type', 'Details', 'Kids_First_Biospecimen_ID', 'sample_id', 'cohort', 'cohort_participant_id'))
    tmpMut <- tmpMut %>%
      group_by(Aberration, Type, Details, variable, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id) %>%
      dplyr::mutate(Variant_Properties = paste0(variable,":", value))  %>%
      group_by(Aberration, Type, Details, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id) %>%
      dplyr::summarise(Variant_Properties = paste(Variant_Properties, collapse = '; ')) %>%
      unique()
  } else {
    tmpMut <- data.frame()
  }
  
  # Add Fusions
  if(!is.null(filtered_fusions)){
    tmpFus <- filtered_fusions %>%
      mutate(Aberration = fusion_name,
             Type = "Fusion",
             Details = paste0("Fusion Type: ", type),
             Variant_Properties = annots) %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
  } else {
    tmpFus <- data.frame()
  }
  
  # Add Copy Number
  if(!is.null(filtered_cnv)){
    tmpCnv <- filtered_cnv %>%
      rowwise() %>%
      mutate(Aberration = hgnc_symbol,
             Type = status,
             Details = paste0("Copy Number: ", copy_number, 
                              " | Pos: ", cytoband),
             Variant_Properties = 'N/A') %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
  } else {
    tmpCnv <- data.frame()
  }
  
  # for RNA data
  if(!is.null(rnaseq_analysis_output)){
    rna_sample <- sample_info %>%
      filter(experimental_strategy == "RNA-Seq") %>%
      dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
  }
  
  # Add Expression
  if(!is.null(rnaseq_analysis_output)){
    tmpExp <- rnaseq_analysis_output$diffexpr.top20 %>%
      rownames_to_column("Aberration") %>%
      mutate(Type = c(rep("Outlier-High (mRNA)", 20), rep("Outlier-Low (mRNA)", 20)),
             Details = paste0("logFC: ", round(logfc, 2)," | TPM: ", tpm),
             Variant_Properties = "N/A",
             cohort  = rna_sample$cohort,
             sample_id = rna_sample$sample_id,
             Kids_First_Biospecimen_ID = rna_sample$Kids_First_Biospecimen_ID,
             cohort_participant_id = rna_sample$cohort_participant_id)  %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
  } else {
    tmpExp <- data.frame()
  }
  
  # Add Pathway
  if(!is.null(rnaseq_analysis_output)){
    tmpPath <- rnaseq_analysis_output$pathways
    tmpPathUp <- tmpPath %>%
      filter(direction == "up") %>%
      mutate(Type = "Pathway Up") %>%
      arrange(padj) %>%
      slice_head(n = 20)
    tmpPathDown <- tmpPath %>%
      filter(direction == "down") %>%
      mutate(Type = "Pathway Down") %>%
      arrange(padj) %>%
      slice_head(n = 20)
    tmpPath <- rbind(tmpPathUp, tmpPathDown)
    tmpPath <- tmpPath %>%
      mutate(Aberration = pathway,
             Details = paste0("Adj. P-Value: ", formatC(padj, format = "e", digits = 2)),
             Variant_Properties = "N/A",
             cohort  = rna_sample$cohort,
             sample_id = rna_sample$sample_id,
             Kids_First_Biospecimen_ID = rna_sample$Kids_First_Biospecimen_ID,
             cohort_participant_id = rna_sample$cohort_participant_id) %>%
      dplyr::select(Aberration, Type, Details, Variant_Properties, Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
  } else {
    tmpPath <- data.frame()
  }
  
  # Merge Together
  allFindingsDF <- unique(rbind(tmpMut, tmpFus, tmpCnv, tmpExp, tmpPath))
  
  # add Ensembl ids and map to targetvalidation.org
  myTable <- allFindingsDF %>%
    left_join(gencode_gtf %>% dplyr::select(gene_id, gene_name), by = c("Aberration" = "gene_name")) %>%
    mutate(TargetValidation = ifelse(is.na(gene_id), "", paste0('<a href = \"https://www.targetvalidation.org/target/',gene_id,'\">',gene_id,"</a>"))) %>%
    dplyr::select(Aberration, Type, Details, Variant_Properties, TargetValidation, 
                  Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id)
  
  return(myTable)
}
