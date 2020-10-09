# get trancriptomically similar patients (top 20 nearest neighbors)

transciptomically_similar <- function(all_cor, clin_data) {
  all_cor <- clin_data %>%
    inner_join(all_cor, by = c("sample_barcode" = "nearest_neighbor")) %>%
    arrange(distance)
  return(all_cor)
}
