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
    filter(Func.refGene %in% c("exonic","splicing","exonic;splicing")) %>% # Func.refGene filter
    mutate(DP = V209, INFO = V219) %>% # unassigned column names
    filter(DP  > 10) %>% # sequencing depth DP > 10
    mutate(AD = gsub("^[0-9][/][0-9]:|:.*", "", INFO)) %>% # capture AD
    separate(AD, c("AD1", "AD2"), sep = ",", extra = "merge", convert = F) %>%
    mutate(AD1 = as.numeric(AD1), 
           AD2 = as.numeric(AD2),
           gnomad211_exome_AF_popmax = as.numeric(gnomad211_exome_AF_popmax),
           gnomad211_genome_AF_popmax = as.numeric(gnomad211_genome_AF_popmax),
           gnomad30_genome_AF = as.numeric(gnomad30_genome_AF)) %>%
    mutate(gnomad211_exome_AF_popmax = replace_na(gnomad211_exome_AF_popmax, 0),
           gnomad211_genome_AF_popmax = replace_na(gnomad211_genome_AF_popmax, 0),
           gnomad30_genome_AF = replace_na(gnomad30_genome_AF, 0)) %>% # replace NA to 0
    filter(AD2 > 0.25 * (AD1 + AD2)) %>% # AD > 0.25 * sum of all ADs
    dplyr::rowwise() %>%
    mutate(AF = max(gnomad211_exome_AF_popmax, gnomad211_genome_AF_popmax, gnomad30_genome_AF)) %>% # calculate AF
    filter(AF < 0.001) %>% # allele frequency < 0.001 in the gnomAD datasets
    mutate(Aberration = AAChange.refGene, 
           Details = paste0("dbSNP: ", V212),
           InterVar_Rank = InterVar_automated) %>% # InterVar_automated as InterVar_Rank
    dplyr::select(Gene.refGene, Aberration, Details, InterVar_Rank, Class, gnomad211_exome_AF_popmax, gnomad211_genome_AF_popmax, gnomad30_genome_AF, DP, INFO) %>%
    unique()
  
  return(germlineOut)
}
