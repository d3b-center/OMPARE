# transcriptome based drug recommendations
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# reference directories
dsigdb_dir <- file.path(ref_dir, "dsigdb")

transcriptome_drug_rec <- function(diffexpr_genes, dsigdb_dat){
  diffexpr_genes <- dsigdb_dat %>%
    inner_join(diffexpr_genes, by = c("Gene" = "genes")) %>%
    rename("Comparison" = "comparison") %>%
    dplyr::select(Drug, Gene, Type, Source, Comparison, logFC)
  return(diffexpr_genes)
}