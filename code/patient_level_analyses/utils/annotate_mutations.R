# annotate mutations vs vus

annotate_mutations <- function(myMutData = mutData, myCancerGenes = cancerGenes) {

  # variant classification filters
  keepVC <- c("Nonsense_Mutation", "Missense_Mutation", 
              "Splice_Region", "Splice_Site",
              "3'UTR", "5'UTR", 
              "In_Frame_Ins", "In_Frame_Del",
              "Frame_Shift_Ins", "Frame_Shift_Del")
  
  # impact filters
  keepVI <- c("MODIFIER", "MODERATE", "HIGH")
  
  # gene to annotation mutation vs vus
  myCancerGenes <- as.character(myCancerGenes$Gene_Symbol)
  
  # filter by biotype, variant class, impact 
  # annotate by cancer gene list (annoFuse)
  mutDataFilt <- v %>%
    filter(BIOTYPE == "protein_coding" &
             Variant_Classification %in% keepVC &
             IMPACT %in% keepVI) %>%
    mutate(Type = ifelse(Hugo_Symbol %in% myCancerGenes, "Mutation", "VUS"))
  
  return(mutDataFilt)
}
