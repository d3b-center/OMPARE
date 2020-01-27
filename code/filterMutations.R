####################################################
# Function to filter mutations-
####################################################

filterMutations <- function(myMutData = mutData, myCancerGenes = cancerGenes) {
  mutDataFilt <- myMutData
  keepVC <- c("Nonsense_Mutation", "Missense_Mutation", "Splice_Region", "3'UTR", "5'UTR", "In_Frame_Del")
  keepVI <- c("MODIFIER", "MODERATE", "HIGH")
  myCancerGenes <- as.character(myCancerGenes$Gene_Symbol)
  
  # Filter by Biotype, Variant class, IMPACT and Cancer gene list
  mutDataFilt <- mutDataFilt %>%
    filter(BIOTYPE == "protein_coding" &
             Variant_Classification %in% keepVC &
             IMPACT %in% keepVI &
             Hugo_Symbol %in% myCancerGenes)
  
  return(mutDataFilt)
}
####################################################
# End Function to filter mutations
####################################################