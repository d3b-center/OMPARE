# get trancriptomically similar patients (top 20 nearest neighbors)

transciptomically_similar <- function(all_cor, clin_data) {
  all_cor <- clin_data %>%
    rownames_to_column("Kids_First_Biospecimen_ID") %>%
    inner_join(all_cor, by = c("Kids_First_Biospecimen_ID" = "nearest_neighbor")) %>%
    arrange(distance)
  return(all_cor)
}
