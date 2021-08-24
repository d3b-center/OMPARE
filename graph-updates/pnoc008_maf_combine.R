
# Downloaded all the pnoc008 maf files - and they are renamed based on subject ID
# **NOTE**:
  # For sample 26, PNOC008-26-7316-8815 is used (not 7316-8812)
  # For sample 5, CHOP sequencing results were used 
  # For sample 2 and 3, WXS results (not WGS) were used

all_maf <- list.files("~/Downloads/maf_pnoc008/")

# read in any vef.maf and any protected.maf to find common colnames
colnames_protected <- readr::read_tsv("~/Downloads/maf_pnoc008/PNOC008-38.consensus_somatic.protected.maf", skip=1) %>%
  colnames()
colnames_protected <- readr::read_tsv("~/Downloads/maf_pnoc008/PNOC008-4.consensus_somatic.vep.maf", skip=1) %>%
  colnames()
common_columns <- intersect(colnames_protected, colnames_protected)

consensus_maf_pnoc008_list <- lapply(all_maf, function(x){
  pnoc008_name = gsub("\\..*", "", x)
  maf_file <- readr::read_tsv(file.path("~/Downloads/maf_pnoc008/",x), skip = 1)
  maf_file <- maf_file %>% 
    select(common_columns) %>%
    mutate(subjectID=pnoc008_name) 
} )

consensus_maf_pnoc008 <- do.call(rbind, consensus_maf_pnoc008_list)%>% distinct()
saveRDS(consensus_maf_pnoc008, "../data/pnoc008_consensus_maf_combined.rds")

################################################ for combining new 

# download the new MAf and rename with the PNOC008 subject ID
new_maf <- list.files("~/Downloads/maf_pnoc008/")
consensus_maf_pnoc008 <- readRDS("../data/pnoc008_consensus_maf_combined.rds")

consensus_maf_new_list <- lapply(new_maf, function(x){
  pnoc008_name = gsub("\\..*", "", x)
  maf_file <- readr::read_tsv(file.path("~/Downloads/maf_pnoc008/",x), skip = 1)
  maf_file <- maf_file %>% 
    mutate(subjectID=pnoc008_name) %>%
    select(colnames(consensus_maf_pnoc008)) 
} )

consensus_maf_new <- do.call(rbind, consensus_maf_new_list)%>% distinct()
consensus_maf_pnoc008 <- rbind(consensus_maf_pnoc008, consensus_maf_new)

saveRDS(consensus_maf_pnoc008, "../data/pnoc008_consensus_maf_combined.rds")

