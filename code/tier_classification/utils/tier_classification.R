# Author: Run Jin
# Classify SNV mutation results into 4 tier [@doi:10.1016/j.jmoldx.2016.10.002]

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

tier_classification <- function(all_findings_output, oncokb_anno, cancer_genes, hotspot_indel, hotspot_snv, cosmic_resistance){
  
  # subset to genes in all_findings_output corresponding to VUS and Mutation
  genes_of_interest <- all_findings_output %>%
    filter(Type %in% c("VUS", "Mutation")) %>% 
    .$Aberration
  hotspot_snv <- hotspot_snv %>% 
    filter(Hugo_Symbol %in% genes_of_interest)
  hotspot_indel <- hotspot_indel %>%
    filter(Hugo_Symbol %in% hotspot_snv$Hugo_Symbol)
  
  # find the oncogenes tumor suppressor genes in OMPARE knowledge base
  oncogenes_tsg <- cancer_genes %>%
    filter(type %in% c("Is.Oncogene", "Oncogene", "Is.Tumor.Suppressor.Gene", "TumorSuppressorGene")) %>%
    .$Gene_Symbol %>% unique
  
  # calculate VAF for subsequent classification 
  oncokb_anno <- oncokb_anno %>% 
    dplyr::mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>% 
    left_join(cosmic_resistance)
  
  # Filter for Tier I
  oncokb_anno_tier1 <- oncokb_anno %>% 
    # filter to variation of interest 
    dplyr::filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', # filter the files 
                                                'Frame_Shift_Del', 'Frame_Shift_Ins', 
                                                'In_Frame_Del', 'In_Frame_Ins', 
                                                'Splice_Region', 'Splice_Site')) %>%
    # filter based on OncoKB annotation level 
    dplyr::filter(HIGHEST_LEVEL %in% c("LEVEL_1","LEVEL_2","LEVEL_3A", "LEVEL_3B", "LEVEL_R1") | HIGHEST_DX_LEVEL %in% c("LEVEL_Dx1","LEVEL_Dx2") | HIGHEST_PX_LEVEL %in% c("LEVEL_Px1","LEVEL_Px2") | cosmic_tier_anno == "Yes" | grepl('\\bpathogenic\\b', CLIN_SIG)) %>%
    # contain only variants that are <1% (less likely to be germline)
    dplyr::filter(gnomAD_AF<0.01| is.na(gnomAD_AF)) %>%
    # filter on VAF over 5% to avoid background noise
    dplyr::filter(VAF>0.05) %>%
    # filter to contain only mutations that has COSMIC ID or has resistance mutation
    dplyr::filter(grepl("COSM", Existing_variation)) %>% 
    # filter to contain only mutations with citations
    dplyr::filter(!is.na(CITATIONS)) %>% 
    # filter to TSG lists in OMPARE database
    dplyr::filter(Hugo_Symbol %in% oncogenes_tsg ) %>% 
    dplyr::mutate(tier_classification = "TIER: I")
  
  
  # Filter for Tier II
  oncokb_anno_tier2 <- oncokb_anno %>% 
    # filter to variation of interest 
    dplyr::filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', # filter the files 
                                                'Frame_Shift_Del', 'Frame_Shift_Ins', 
                                                'In_Frame_Del', 'In_Frame_Ins', 
                                                'Splice_Region', 'Splice_Site')) %>%
    # filter based on OncoKB annotation level 
    dplyr::filter(HIGHEST_LEVEL %in% c("LEVEL_4","LEVEL_R2") | HIGHEST_DX_LEVEL %in% c("LEVEL_Dx3") | HIGHEST_PX_LEVEL %in% c("LEVEL_Px3") | cosmic_tier_anno == "No" | (grepl('likely_pathogenic', CLIN_SIG) & !grepl('\\bpathogenic\\b', CLIN_SIG))) %>%
    # contain only variants that are <1% (less likely to be germline)
    dplyr::filter(gnomAD_AF<0.01| is.na(gnomAD_AF)) %>%
    # filter on VAF over 5% to avoid background noise
    dplyr::filter(VAF>0.05) %>%
    # filter to contain only mutations that has COSMIC ID or has resistance mutation
    dplyr::filter(grepl("COSM", Existing_variation)) %>% 
    # filter to contain only mutations with citations
    dplyr::filter(!is.na(CITATIONS)) %>% 
    # filter to TSG lists in OMPARE database
    dplyr::filter(Hugo_Symbol %in% oncogenes_tsg ) %>% 
    dplyr::mutate(tier_classification = "TIER: II")
  
  
  # Filter for Tier III
  oncokb_anno_tier3 <- oncokb_anno %>% 
    # filter to variation of interest 
    dplyr::filter(Variant_Classification %in% c('Missense_Mutation',  # filter the files 
                                                'In_Frame_Del', 'In_Frame_Ins')) %>%
    # filter to samples that do not have OncoKB annotation 
    dplyr::filter(is.na(HIGHEST_LEVEL)) %>%
    dplyr::filter(is.na(HIGHEST_DX_LEVEL)) %>%
    dplyr::filter(is.na(HIGHEST_PX_LEVEL)) %>%
    # contain only variants that are <1% (less likely to be germline)
    dplyr::filter(gnomAD_AF<0.01| is.na(gnomAD_AF)) %>%
    # filter on VAF over 5% to avoid background noise
    dplyr::filter(VAF>0.05) %>%
    dplyr::mutate(tier_classification = "TIER: III")
  
  # Filter for Tier IV
  oncokb_anno_tier4 <- oncokb_anno %>% 
    # filter to variation of interest 
    dplyr::filter(Variant_Classification =='Missense_Mutation') %>%
    # filter to samples that do not have OncoKB annotation or related citations
    dplyr::filter(is.na(HIGHEST_LEVEL)) %>%
    dplyr::filter(is.na(HIGHEST_DX_LEVEL)) %>%
    dplyr::filter(is.na(HIGHEST_PX_LEVEL)) %>%
    dplyr::filter(is.na(CITATIONS)) %>% 
    # contain only variants that are <1% (less likely to be germline)
    dplyr::filter(gnomAD_AF>=0.01) %>%
    # filter on VAF over 5% to avoid background noise
    dplyr::filter((VAF>=0.4 & VAF<=0.6) | VAF>=0.9) %>%
    dplyr::filter(!grepl("deleterious", SIFT)) %>%
    dplyr::filter(!grepl("damaging", PolyPhen)) %>%
    dplyr::mutate(tier_classification = "TIER: IV")
  
  # combine all the variants in MAF that contain Tier information 
  oncokb_anno_tiers <- bind_rows(oncokb_anno_tier1, oncokb_anno_tier2, oncokb_anno_tier3, oncokb_anno_tier4) %>%
    dplyr::select(Hugo_Symbol, Variant_Classification, HGVSp_Short, tier_classification)
  
  
  ######################### annotate MSK SNV hotspot database to all findings table
  
  # annotate SNV first
  if(nrow(hotspot_snv)>0){
    # modify the clinical output files for easier manipulation
    all_findings_output_snv <- all_findings_output %>% 
      dplyr::filter(grepl("Missense_Mutation", Details) | grepl("Nonsense_Mutation", Details) | grepl("Splice_Site", Details) | grepl("Splice_Region", Details)) %>%
      dplyr::mutate(position = gsub(".*HGVSp: p..", "", Details)) %>% 
      dplyr::mutate(position =str_sub(position, end = -2))
    
    all_findings_output_hotspots <- data.frame()
    
    snv_location <- hotspot_snv %>% 
      dplyr::select(Hugo_Symbol, Amino_Acid_Position) %>% 
      distinct()
    
    for(p in 1:nrow(snv_location)){
      
      each_entry <- snv_location[p,]
      hotspot_gene <- each_entry$Hugo_Symbol
      hotspot_aa_position <- each_entry$Amino_Acid_Position
      
      # find all the specific AA mutations in the hotspot table
      hotspot_aa <- hotspot_snv %>% 
        dplyr::filter(Amino_Acid_Position == hotspot_aa_position,
                      Hugo_Symbol==hotspot_gene) %>%
        dplyr::pull(HGVSp_short)
      
      # reformat to make it matchable 
      hotspot_aa<-paste(noquote(hotspot_aa),collapse='|')
      
      # first deal with all clinical findings table
      all_findings_output_hotspot <- all_findings_output_snv %>%
        dplyr::filter(Aberration == hotspot_gene,
                      position == hotspot_aa_position) %>%
        dplyr::select(-position)
      
      if(nrow(all_findings_output_hotspot)>0){
        all_findings_output_hotspot <- all_findings_output_hotspot %>%
          dplyr::mutate(Variant_Properties= case_when(
            grepl(hotspot_aa, Details) ~ paste0("Cancer Hotspot; ", Variant_Properties),
            TRUE ~ paste0("Cancer Hotspot Location; ", Variant_Properties)
          ))
      }
      all_findings_output_hotspots <- bind_rows(all_findings_output_hotspots, all_findings_output_hotspot)
    }
    
    # add the annotations to the file 
    all_findings_output_not_hotspot <- all_findings_output %>% 
      dplyr::filter(!Details %in% all_findings_output_hotspots$Details)
    all_findings_output <- bind_rows(all_findings_output_not_hotspot, all_findings_output_hotspots)
  }
  
  ######################### annotate MSK Indel hotspot database to all findings table
  
  # annotate Indel now
  if(nrow(hotspot_indel)>0){
    all_findings_output_indels <- data.frame()
    
    all_findings_output_mnv <- all_findings_output %>%
      dplyr::filter(grepl("Frame_Shift_Del", Details) | grepl("Frame_Shift_Ins", Details) | grepl("In_Frame_Del", Details) | grepl("In_Frame_Ins", Details)) %>%
      dplyr::mutate(specific = gsub(".*HGVSp: ", "", Details)) 
    
    for(q in 1:nrow(hotspot_indel)){
      
      each_indel <- hotspot_indel[q,]
      hotspot_gene_indel <- each_indel$Hugo_Symbol
      hotspot_specific_indel <- each_indel$HGVSp_short
      
      all_findings_output_indel <- all_findings_output_mnv %>%
        dplyr::filter(Aberration == hotspot_gene_indel) %>%
        dplyr::filter(specific == hotspot_specific_indel) %>%
        dplyr::mutate(Variant_Properties= paste0("Cancer Hotspot; ", Variant_Properties)) %>%
        dplyr::select(-specific)
      all_findings_output_indels <- bind_rows(all_findings_output_indel, all_findings_output_indels)
      
    }
    
    # add the annotations to the file 
    all_findings_output_no_indels <- all_findings_output %>% 
      dplyr::filter(!Details %in% all_findings_output_indels$Details)
    all_findings_output <- bind_rows(all_findings_output_no_indels, all_findings_output_indels)
  }
  
  
  ######################### annotate TIER information to all findings output
  
  # annotate TIER info
  if(nrow(oncokb_anno_tiers)){
    all_findings_output_tiers <- data.frame()
    
    for(i in 1:nrow(oncokb_anno_tiers)) {
      # for each mutation, find the matching entries in the all findings table (if present)
      each <- oncokb_anno_tiers[i,]
      gene_symbol <- each$Hugo_Symbol
      variant_classification <- each$Variant_Classification
      hgvs_short <- each$HGVSp_Short
      tier_class <- each$tier_classification
      
      # make changes to variant property column when TIER information is available
      all_findings_output_tier <- all_findings_output %>% 
        dplyr::filter(Aberration == gene_symbol & grepl(variant_classification, Details) & grepl(hgvs_short, Details)) %>%
        dplyr::mutate(Variant_Properties =paste0(tier_class, "; ", Variant_Properties))
      
      all_findings_output_tiers <- rbind(all_findings_output_tiers, all_findings_output_tier)
    }
    
    # filter to the ones that do not have modification
    all_findings_output_no_tier <- all_findings_output %>% 
      dplyr::filter(!Details %in% all_findings_output_tiers$Details)
    all_findings_output <- bind_rows(all_findings_output_no_tier, all_findings_output_tiers)
  }
  
  # return annotated output
  return(all_findings_output)
}
