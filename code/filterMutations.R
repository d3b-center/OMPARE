####################################################
# Function to filter mutations-
####################################################
filterMutations <- function(myMutData=mutData, myCancerGenes=cancerGenes) {
  mutDataFilt <- myMutData
  
  #Filter to only protein coding genes
  mutDataFilt <- mutDataFilt[mutDataFilt[,"BIOTYPE"]=="protein_coding",]
  
  #Filter by Variant Classification
  keepVC <- c("Missense_Mutation", "Splice_Region", "3'UTR", "5'UTR", "In_Frame_Del")
  mutDataFilt <- mutDataFilt[mutDataFilt[,"Variant_Classification"]%in%keepVC,]
  
  #Filter by Variant IMPACT
  keepVI <- c("MODIFIER", "MODERATE", "HIGH")
  mutDataFilt <- mutDataFilt[mutDataFilt[,"IMPACT"]%in%c(keepVI),]
  
  #Filter by Cancer Gene Census
  myCancerGenes <- as.character(myCancerGenes[,1])
  mutDataFilt <- mutDataFilt[mutDataFilt[,"Hugo_Symbol"]%in%myCancerGenes,]
  return(mutDataFilt)
}
####################################################
# End Function to filter mutations
####################################################