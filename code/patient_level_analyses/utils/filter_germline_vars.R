# code to filter germline variants
filter_germline_vars <- function(mutData.germ = mutData.germ, myGermlineMarkers = germline_markers){
  
  # group by genes
  myGermlineMarkers <- myGermlineMarkers %>% 
    group_by(Gene) %>%
    summarise(Class = toString(Class))
  
  germlineData <- mutData.germ
  
  # Filter only to germlineMarkers genes
  germlineData <- merge(germlineData, myGermlineMarkers, by.x = 'Gene.refGene', by.y = 'Gene')
  
  # Additional Filters
  germlineOut  <- germlineData %>% 
    # Func.refGene filter
    filter(Func.refGene %in% c("exonic","splicing","exonic;splicing")) %>% 
    # unassigned column names
    mutate(DP = V209, INFO = V219) %>%
    # sequencing depth DP > 10
    filter(DP  > 10) %>% 
    # capture AD
    mutate(AD = gsub("^[0-9][/][0-9]:|:.*", "", INFO)) %>% 
    separate(AD, c("AD1", "AD2"), sep = ",", extra = "merge", convert = F) %>%
    mutate(AD1 = as.numeric(AD1), 
           AD2 = as.numeric(AD2),
           gnomad211_exome_AF_popmax = as.numeric(gnomad211_exome_AF_popmax),
           gnomad211_genome_AF_popmax = as.numeric(gnomad211_genome_AF_popmax),
           gnomad30_genome_AF = as.numeric(gnomad30_genome_AF)) %>%
    mutate(gnomad211_exome_AF_popmax = replace_na(gnomad211_exome_AF_popmax, 0),
           gnomad211_genome_AF_popmax = replace_na(gnomad211_genome_AF_popmax, 0),
           # replace NA to 0
           gnomad30_genome_AF = replace_na(gnomad30_genome_AF, 0)) %>% 
    # AD > 0.25 * sum of all ADs
    filter(AD2 > 0.25 * (AD1 + AD2)) %>% 
    dplyr::rowwise() %>%
    # calculate AF
    mutate(AF = max(gnomad211_exome_AF_popmax, gnomad211_genome_AF_popmax, gnomad30_genome_AF)) %>% 
    # allele frequency < 0.001 in the gnomAD datasets
    filter(AF < 0.001) %>% 
    mutate(Aberration = AAChange.refGene, 
           Details = paste0("dbSNP: ", V212),
           # add Clinvar prediction
           Clinvar_Pred = CLNSIG, 
           # InterVar_automated as InterVar_Rank
           InterVar_Rank = InterVar_automated) %>% 
    dplyr::select(Gene.refGene, Aberration, Details, Clinvar_Pred, InterVar_Rank, Class, gnomad211_exome_AF_popmax, gnomad211_genome_AF_popmax, gnomad30_genome_AF, DP, INFO) %>%
    unique() %>%
    as.data.frame()
  
  return(germlineOut)
}
