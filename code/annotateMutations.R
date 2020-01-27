####################################################
# Function to annotate VUS mutations
####################################################

annotateMutations <- function(myMutData = mutData, myCancerGenes = cancerGenes) {
  mutDataFilt <- myMutData
  keepVC <- c("Missense_Mutation", "Splice_Region", "3'UTR", "5'UTR", "In_Frame_Del")
  keepVI <- c("MODIFIER", "MODERATE", "HIGH")
  myCancerGenes <- as.character(myCancerGenes$Gene_Symbol)
  
  # Filter by Biotype, Variant class, IMPACT 
  # Annotate by Cancer gene list
  mutDataFilt <- mutDataFilt %>%
    filter(BIOTYPE == "protein_coding" &
             Variant_Classification %in% keepVC &
             IMPACT %in% keepVI) %>%
    mutate(Type = ifelse(Hugo_Symbol %in% myCancerGenes, "Mutation", "VUS"))
  
  return(mutDataFilt)
}
####################################################
# End Function to annotate VUS mutations
####################################################