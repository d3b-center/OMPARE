############################
# Purpose: Germline Analysis
############################

germlineAnalysis <- function(mutData = mutData.germ, myGermlineMarkers = germlineMarkers){
  
  # group by genes
  myGermlineMarkers <- myGermlineMarkers %>% 
    group_by(Gene) %>%
    summarise(Class = toString(Class))
  
  germlineData <- mutData

  # Filter only to germlineMarkers genes
  germlineData <- merge(germlineData, myGermlineMarkers, by.x = 'Gene.refGene', by.y = 'Gene')
  
  # Additional Filters
  germlineOut  <- germlineData %>% 
    mutate(DP = V209, INFO = V219) %>%
    filter(DP  > 10) %>% # sequencing depth DP > 10
    filter(AAChange.refGene != ".") %>% # AAChange.refGene not == "."
    mutate(AD = gsub("^[0-9][/][0-9]:|:.*", "", INFO)) %>%
    separate(AD, c("AD1", "AD2"), sep = ",", extra = "merge", convert = F) %>%
    mutate(AD1 = as.numeric(AD1), 
           AD2 = as.numeric(AD2),
           gnomad211_exome_AF = as.numeric(gnomad211_exome_AF),
           gnomad211_genome_AF = as.numeric(gnomad211_genome_AF),
           gnomad30_genome_AF = as.numeric(gnomad30_genome_AF)) %>%
    filter(AD2 > 0.25 * (AD1 + AD2)) %>% # AD > 0.25 * sum of all AD's
    dplyr::rowwise() %>%
    mutate(AF = min(gnomad211_exome_AF, gnomad211_genome_AF, gnomad30_genome_AF,  na.rm = TRUE)) %>% # calculate AF
    filter(AF < 0.01) %>% # allele frequency < 0.01 in the gnomAD datasets
    mutate(Aberration = AAChange.refGene, 
           Details = paste0("dbSNP: ", V212),
           InterVar_Rank = InterVar_automated) %>%
    dplyr::select(Gene.refGene, Aberration, Details, InterVar_Rank, Class, gnomad211_exome_AF, gnomad211_genome_AF, gnomad30_genome_AF, DP, INFO) %>%
    unique()
  
  # Filter to only genes in germline markers
  # germlineMarkersGenes <- unique(germlineMarkers[,"Gene"]);
  # germlineData <- germlineData[germlineData[,1] %in% germlineMarkersGenes,]
  
  # Filter by annotation
  # germlineData <- germlineData[germlineData[,"Consequence"]!="synonymous_variant",]
  # germlineData <- germlineData[grepl("pathogenic", germlineData[,"CLIN_SIG"]),]
  # germlineOut <- merge(germlineMarkers, germlineData, by.y ="Hugo_Symbol", by.x="Gene")
  # germlineOut <- germlineOut[,c("Disease", "Gene", "HGVSp_Short")]
  return(germlineOut)
}
