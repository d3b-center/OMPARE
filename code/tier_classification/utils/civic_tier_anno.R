# Author: Run Jin
# Annotate with CIVIC information 
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(data.table)
})

civic_annotation <- function(all_findings_output, civic_ref_dir, civic_output){
  # give each entry an ID for easier match
  all_findings_output <- all_findings_output %>%
    tibble::rownames_to_column("findings_id")
  
  # subset for different type
  all_findings_output_mut <- all_findings_output %>%
    dplyr::filter(Type %in% c("VUS", "Mutation"))
  all_findings_output_cnv <- all_findings_output %>% 
    dplyr::filter(Type %in% c("Gain", "Amplification", "Loss"))
  all_findings_output_fusion <- all_findings_output %>%
    dplyr::filter(Type == "Fusion")
  all_findings_output_expression <- all_findings_output %>%
    dplyr::filter(Type %in% c("Outlier-High (mRNA)",
                              "Outlier-Low (mRNA)"))
  
  ########## First deal with specific mutation 
  c_dot_civic <- read_tsv(file.path(civic_ref_dir, "c_dot_civic.tsv")) %>%
    dplyr::filter(grepl("\\(", variant)) %>%
    dplyr::mutate(variant = gsub(" \\s*\\([^\\)]+\\)", "", variant)) 
  snv_remaining <- read_tsv(file.path(civic_ref_dir, "snv_remaining_civic.tsv"),
                            guess_max = 10000) %>% 
    dplyr::mutate(gene = case_when(
      gene == "H3-3A" ~ "H3F3A",
      TRUE ~ gene
    ))
  
  # bind rows and reformat
  snv_combined <- rbind(c_dot_civic, snv_remaining) %>%
    dplyr::mutate(variant = gsub("FS", "fs", variant),
                  variant = gsub("Ter", "\\*", variant),
                  variant = gsub("DUP", "dup", variant),
                  variant = gsub("DEL", "del", variant),
                  variant = gsub("INS", "ins", variant))
  
  ######################### annotate specific mutations -------------------------------
  # define df for taking the hit 
  snv_specific_hit <- data.frame()
  for(i in 1:nrow(snv_combined)){
    gene_of_interest <- snv_combined[i,2]
    variant_of_interest <- snv_combined[i,4]
    
    hit_each <- all_findings_output_mut %>%
      dplyr::filter(Aberration == gene_of_interest,
                    grepl(variant_of_interest, "Details")) %>%
      dplyr::left_join(snv_combined[i,], by = c("Aberration" = "gene"))
    snv_specific_hit <- bind_rows(snv_specific_hit, hit_each)
  }
  # find the un-annotated one
  all_findings_output_mut_remaining <- all_findings_output_mut %>%
    dplyr::filter(!findings_id %in% snv_specific_hit$findings_id)
  
  ######################### annotate exon and domain mutations -------------------------------
  # first find the domain and exon genes that are not in the general dataset already 
  mutation_exon_civic <- read_tsv(file.path(civic_ref_dir, "mutation_exon_civic_annotated.tsv"),
                                  guess_max = 10000) 
  
  # annotate the domain information
  mutation_domain_civic <- read_tsv(file.path(civic_ref_dir, "mutation_domain_civic_annotated.tsv"),
                                    guess_max = 10000) 
  
  mutation_domain_exon_civic <- bind_rows(mutation_exon_civic, 
                                          mutation_domain_civic)
  # define dataframe to store results
  rm(hit_each)
  snv_exon_domain_hit <- data.frame()
  for(i in 1:nrow(mutation_domain_exon_civic)){
    # find the gene, and start and end 
    gene_of_interest <- mutation_domain_exon_civic[i,1]
    protein_start <- mutation_domain_exon_civic[1,3]
    protein_end <- mutation_domain_exon_civic[1,4]
    
    # find each hit
    hit_each <- all_findings_output_mut_remaining %>%
      dplyr::filter(Aberration == gene_of_interest) %>%
      dplyr::mutate(aa_location = gsub("fs*.*", "", Details)) %>%
      dplyr::mutate(aa_location = as.numeric(str_extract(Details, "[0-9]+"))) %>%
      dplyr::filter(aa_location >= protein_start & 
                      aa_location <= protein_end) %>%
      dplyr::left_join(mutation_domain_exon_civic[i,], by = c("Aberration" = "gene"))
    
    # combine all the hits together
    snv_exon_domain_hit <- bind_rows(snv_exon_domain_hit, hit_each)
  }
  snv_exon_domain_hit <- snv_exon_domain_hit %>%
    dplyr::select(-c("aa_location", "protein_start", "protein_end"))
  # find the un-annotated one
  all_findings_output_mut_remaining <- all_findings_output_mut_remaining %>%
    dplyr::filter(!findings_id %in% snv_exon_domain_hit$findings_id)
  
  ######################### annotate exon insertion and deletions mutations -------------------------------
  deletion_exon <- read_tsv(file.path(civic_ref_dir, "del_exon_civic_annotated.tsv"))
  insertion_exon <- read_tsv(file.path(civic_ref_dir, "ins_exon_civic_annotated.tsv"))
  
  # define dataframe to store results
  rm(hit_each)
  del_exon_hit <- data.frame()
  for(i in 1:nrow(deletion_exon)){
    # find the gene, and start and end 
    gene_of_interest <- deletion_exon[i,1]
    protein_start <- deletion_exon[1,3]
    protein_end <- deletion_exon[1,4]
    
    # find each hit
    hit_each <- all_findings_output_mut_remaining %>%
      dplyr::filter(Aberration == gene_of_interest) %>%
      dplyr::filter(grepl("_Del", Details)) %>% 
      dplyr::mutate(aa_location = gsub("fs*.*", "", Details)) %>%
      dplyr::mutate(aa_location = as.numeric(str_extract(Details, "[0-9]+"))) %>%
      dplyr::filter(aa_location >= protein_start & 
                      aa_location <= protein_end) %>%
      dplyr::left_join(deletion_exon[i,], by = c("Aberration" = "gene"))
    
    # combine all the hits together
    del_exon_hit <- bind_rows(del_exon_hit, hit_each)
  }
  del_exon_hit <- del_exon_hit %>%
    dplyr::select(-c("aa_location", "protein_start", "protein_end"))
  
  # find the un-annotated one
  all_findings_output_mut_remaining <- all_findings_output_mut_remaining %>%
    dplyr::filter(!findings_id %in% del_exon_hit$findings_id)
  
  # define dataframe to store results
  rm(hit_each)
  ins_exon_hit <- data.frame()
  for(i in 1:nrow(insertion_exon)){
    # find the gene, and start and end 
    gene_of_interest <- insertion_exon[i,1]
    protein_start <- insertion_exon[1,3]
    protein_end <- insertion_exon[1,4]
    
    # find each hit
    hit_each <- all_findings_output_mut_remaining %>%
      dplyr::filter(Aberration == gene_of_interest) %>%
      dplyr::filter(grepl("_Ins", Details)) %>% 
      dplyr::mutate(aa_location = gsub("fs*.*", "", Details)) %>%
      dplyr::mutate(aa_location = as.numeric(str_extract(Details, "[0-9]+"))) %>%
      dplyr::filter(aa_location >= protein_start & 
                      aa_location <= protein_end) %>%
      dplyr::left_join(insertion_exon[i,], by = c("Aberration" = "gene"))
    
    # combine all the hits together
    ins_exon_hit <- bind_rows(ins_exon_hit, hit_each)
  }
  ins_exon_hit <- ins_exon_hit %>%
    dplyr::select(-c("aa_location", "protein_start", "protein_end"))
  
  # find the un-annotated one
  all_findings_output_mut_remaining <- all_findings_output_mut_remaining %>%
    dplyr::filter(!findings_id %in% ins_exon_hit$findings_id)
  
  ######################### annotate general mutation and combine results -------------------------------
  # finally annotate to general hits 
  mutation_general <- read_tsv(file.path(civic_ref_dir, "mutation_general_civic.tsv"),
                               guess_max = 10000) %>% 
    dplyr::mutate(gene = case_when(
      gene == "H3-3A" ~ "H3F3A",
      TRUE ~ gene
    ))
  mutation_general_hit <- all_findings_output_mut_remaining %>% 
    dplyr::filter(Aberration %in% mutation_general$gene) %>%
    dplyr::left_join(mutation_general, by = c("Aberration" = "gene"))
  
  # combine all mutation hits 
  mutation_hit <- rbind(snv_specific_hit,
                        snv_exon_domain_hit,
                        del_exon_hit,
                        ins_exon_hit,
                        mutation_general_hit) %>%
    as.data.frame()
  
  # find the un-annotated one
  mutation_non_hit <- all_findings_output_mut %>%
    dplyr::filter(!findings_id %in% mutation_hit$findings_id)
  
  ################# Annotate expression -------------------------------
  expression_all_civic <- read_tsv(file.path(civic_ref_dir, "expression_all_civic.tsv"))
  expression_hit <- all_findings_output_expression %>%
    dplyr::filter(Aberration %in% expression_all_civic$gene) %>%
    dplyr::left_join(expression_all_civic, by = c("Aberration" = "gene")) %>%
    as.data.frame()
  
  expression_non_hit <- all_findings_output_expression %>%
    dplyr::filter(!Aberration %in% expression_all_civic$gene)
  
  ################ Annotate fusion ----------------------------------
  # read in necessary file
  fusion_column_only <- read_tsv(file.path(civic_ref_dir, "fusion_column_only.tsv")) %>% 
    dplyr::mutate(variant_match = gsub("::", "_", variant))
  fusion_column_hit <- all_findings_output_fusion %>%
    dplyr::filter(Aberration %in% fusion_column_only$variant_match) %>%
    dplyr::left_join(fusion_column_only, by = c("Aberration" = "gene")) %>%
    dplyr::select(-c("variant_match")) %>%
    as.data.frame()
  
  fusion_word_only <- read_tsv(file.path(civic_ref_dir, "fusion_word_only.tsv"))
  all_findings_output_fusion_sep <- all_findings_output_fusion %>%
    dplyr::filter(!findings_id %in% fusion_column_hit$variant_id) %>%
    tidyr::separate(Aberration, c("geneA", "geneB")) 
  
  # merge with CIVIC information based on gene 
  fusion_word_geneA <- all_findings_output_fusion_sep %>%
    dplyr::filter(geneA %in% fusion_word_only$gene) %>%
    dplyr::left_join(fusion_word_only, by = c("geneA" = "gene")) 
  fusion_word_geneB <- all_findings_output_fusion_sep %>%
    dplyr::filter(geneB %in% fusion_word_only$gene) %>%
    dplyr::left_join(fusion_word_only, by = c("geneB" = "gene")) 
  
  # combine together to get the hit
  fusion_word_hit <- bind_rows(fusion_word_geneA, fusion_word_geneB) %>%
    dplyr::mutate(Aberration = paste0(geneA, "_", geneB)) %>%
    dplyr::select(-c("geneA", "geneB")) %>%
    as.data.frame()
  
  # get all fusions 
  fusion_hit <- rbind(fusion_word_hit, fusion_column_hit)
  fusion_non_hit <- all_findings_output_fusion %>%
    dplyr::filter(!findings_id %in% fusion_hit$findings_id)
  
  ################ Annotate amplification and loss ----------------------------------
  # read in necessary file
  amplification_only <- read_tsv(file.path(civic_ref_dir, "amplification_only_civic.tsv"))
  cnv_loss <- read_tsv(file.path(civic_ref_dir, "cnv_loss_civic.tsv"))
  deletion_general_civic <- read_tsv(file.path(civic_ref_dir, "deletion_general_civic.tsv"))
  combined_loss <- rbind(cnv_loss, deletion_general_civic)
  
  # annotate amplification first
  amplification_hit <- all_findings_output_cnv %>%
    dplyr::filter(Aberration %in% amplification_only$gene) %>%
    dplyr::filter(Type %in% c("Gain", "Amplification")) %>%
    dplyr::left_join(amplification_only, by = c("Aberration" = "gene")) %>%
    as.data.frame()
  
  # then annotate loss
  loss_hit <- all_findings_output_cnv %>%
    dplyr::filter(Aberration %in% combined_loss$gene) %>%
    dplyr::filter(Type %in% c("Loss")) %>%
    dplyr::left_join(combined_loss, by = c("Aberration" = "gene")) %>%
    as.data.frame()
  
  # find the hit and non hit
  cnv_hit <- rbind(amplification_hit, loss_hit)
  cnv_non_hit <- all_findings_output_cnv %>%
    dplyr::filter(!findings_id %in% cnv_hit$findings_id)
  
  ################ Annotate alteration ----------------------------------
  # combine all non hit
  all_non_hit <- rbind(mutation_non_hit,
                       expression_non_hit,
                       fusion_non_hit,
                       cnv_non_hit) 
  # read in alteration
  alteration_civic <- read_tsv(file.path(civic_ref_dir, "alteration_civic.tsv"),
                               guess_max = 10000) 
  
  # from the non-hit, rescue all alteration
  alteration_hit <- all_non_hit %>%
    dplyr::filter(Aberration %in% alteration_civic$gene) %>%
    as.data.frame()
  
  # combined all hit together
  all_findings_annotated <- rbind(mutation_hit,
                                  expression_hit,
                                  fusion_hit,
                                  cnv_hit,
                                  alteration_hit)
  all_findings_non_annotated <- all_findings_output %>%
    dplyr::filter(!findings_id %in% all_findings_annotated$findings_id)
  
  # add them together and write out
  all_findings_civic_annotated <- bind_rows(all_findings_annotated, 
                                            all_findings_non_annotated)
  all_findings_civic_annotated %>% 
    write_tsv(file.path(civic_output, paste0("all_findings_output_civic_annotated.tsv")))
  
  # now find new files
  all_findings_civic_annotated_clean <- all_findings_civic_annotated %>% 
    dplyr::select("findings_id", "Aberration", "Type", "Details", "Variant_Properties", "TargetValidation", "evidence_level",
                  "Kids_First_Biospecimen_ID", "sample_id", "cohort", "cohort_participant_id") %>%
    dplyr::group_by(Details) %>%
    dplyr::mutate(evidence_level_civic = paste0(unique(evidence_level), collapse = ",")) %>%
    dplyr::mutate(evidence_level_civic = gsub("NA,", "", evidence_level_civic)) %>%
    dplyr::mutate(evidence_level_civic = gsub("NA", "", evidence_level_civic)) %>%
    dplyr::select(-c("evidence_level")) %>% 
    distinct() %>%
    dplyr::mutate(Variant_Properties = case_when(
      grepl("A", evidence_level_civic) | grepl("B", evidence_level_civic) ~ paste0("TIER: I;", Variant_Properties),
      grepl("C", evidence_level_civic) | grepl("D", evidence_level_civic) ~ paste0("TIER: II;", Variant_Properties),
      grepl("E", evidence_level_civic) ~ paste0("TIER: III;", Variant_Properties),
      TRUE ~ Variant_Properties
    )) %>%
    # remove potential differences
    dplyr::mutate(Variant_Properties = gsub("TIER: I;TIER: I;", "TIER: I;", Variant_Properties),
                  Variant_Properties = gsub("TIER: I;TIER: II;", "TIER: I;", Variant_Properties),
                  Variant_Properties = gsub("TIER: I;TIER: III;", "TIER: I;", Variant_Properties),
                  Variant_Properties = gsub("TIER: I;TIER: IV;", "TIER: I;", Variant_Properties),
                  Variant_Properties = gsub("TIER: II;TIER: II;", "TIER: II;", Variant_Properties),
                  Variant_Properties = gsub("TIER: II;TIER: III;", "TIER: II;", Variant_Properties),
                  Variant_Properties = gsub("TIER: II;TIER: IV;", "TIER: II;", Variant_Properties),
                  Variant_Properties = gsub("TIER: III;TIER: III;", "TIER: III;", Variant_Properties),
                  Variant_Properties = gsub("TIER: IV;TIER: IV;", "TIER: IV;", Variant_Properties),
                  Variant_Properties = gsub("TIER: III;TIER: IV;", "TIER: III;", Variant_Properties),
                  # also select when flipped
                  Variant_Properties = gsub("TIER: II;TIER: I;", "TIER: I;", Variant_Properties),
                  Variant_Properties = gsub("TIER: III;TIER: I;", "TIER: I;", Variant_Properties),
                  Variant_Properties = gsub("TIER: IV;TIER: I;", "TIER: I;", Variant_Properties),
                  Variant_Properties = gsub("TIER: III;TIER: II;", "TIER: II;", Variant_Properties),
                  Variant_Properties = gsub("TIER: IV;TIER: II;", "TIER: II;", Variant_Properties),
                  Variant_Properties = gsub("TIER: IV;TIER: III;", "TIER: III;", Variant_Properties)) %>%
    dplyr::mutate(Variant_Properties = gsub(";NA", "", Variant_Properties)) %>%
    dplyr::select(-c("findings_id", "evidence_level_civic"))
  
  # write out file
  return(all_findings_civic_annotated_clean)
}

