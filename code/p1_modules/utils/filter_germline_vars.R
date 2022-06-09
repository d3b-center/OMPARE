# code to filter germline variants
filter_germline_vars <- function(patient_dir = patient_dir, germline_markers = germline_markers){
  
  # group by genes
  germline_markers <- germline_markers %>% 
    group_by(Gene) %>%
    dplyr::summarise(Class = toString(Class))
  
  # input germline data (this is a big file so no need to read every time)
  mutFiles <- list.files(path = file.path(patient_dir, 'simple-variants'), pattern = 'hg38_multianno.txt.gz', recursive = TRUE, full.names = T)
  germlineData <- data.table::fread(mutFiles, stringsAsFactors = F)
  
  # Filter only to germlineMarkers genes
  germlineData <- germlineData %>%
    inner_join(germline_markers, by = c("Gene.refGene" = "Gene"))
  
  # Additional Filters
  germlineOut  <- germlineData %>% 
    # Func.refGene filter
    filter(Func.refGene %in% c("exonic","splicing","exonic;splicing")) %>% 
    # unassigned column names
    dplyr::mutate(DP = V209, INFO = V219) %>%
    # sequencing depth DP > 10
    filter(DP  > 10) %>% 
    # capture AD
    dplyr::mutate(AD = gsub("^[0-9][/][0-9]:|:.*", "", INFO)) %>% 
    separate(AD, c("AD1", "AD2"), sep = ",", extra = "merge", convert = F) %>%
    dplyr::mutate(AD1 = as.numeric(AD1), 
                  AD2 = as.numeric(AD2),
                  gnomad211_exome_AF_popmax = as.numeric(gnomad211_exome_AF_popmax),
                  gnomad211_genome_AF_popmax = as.numeric(gnomad211_genome_AF_popmax),
                  gnomad30_genome_AF = as.numeric(gnomad30_genome_AF)) %>%
    dplyr::mutate(gnomad211_exome_AF_popmax = replace_na(gnomad211_exome_AF_popmax, 0),
                  gnomad211_genome_AF_popmax = replace_na(gnomad211_genome_AF_popmax, 0),
                  # replace NA to 0
                  gnomad30_genome_AF = replace_na(gnomad30_genome_AF, 0)) %>% 
    # AD > 0.25 * sum of all ADs
    filter(AD2 > 0.25 * (AD1 + AD2)) %>% 
    dplyr::rowwise() %>%
    # calculate AF
    dplyr::mutate(AF = max(gnomad211_exome_AF_popmax, gnomad211_genome_AF_popmax, gnomad30_genome_AF)) %>% 
    # allele frequency < 0.001 in the gnomAD datasets
    filter(AF < 0.001) %>% 
    dplyr::mutate(Aberration = AAChange.refGene, 
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