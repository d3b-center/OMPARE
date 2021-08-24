# transcriptome based drug recommendations
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

transcriptome_drug_rec <- function(diffexpr_genes){
  
  # get gene symbols
  genes <- unique(diffexpr_genes$genes)
  
  # obtained from ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/chembl_29_sqlite.tar.gz
  chembldb <- file.path(ref_dir, "chembl", "chembl_29_sqlite", "chembl_29.db")
  resultsPath <- system.file("extdata", "results", package = "drugTargetInteractions")
  config <- genConfig(chemblDbPath = chembldb, resultsPath = resultsPath)
  downloadUniChem(config = config)
  cmpIdMapping(config = config)
  
  # generate a mapping of uniprot ids and gene symbols only once save as rds
  map_file <- file.path(ref_dir, "uniprot_genesymbol_map.rds")
  if(!file.exists(map_file)){
    hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
    annotation.uni <- getBM(attributes=c("uniprotswissprot", "hgnc_symbol"), mart=hsmart)
    saveRDS(annotation.uni, file = map_file)
  } else {
    annotation.uni <- readRDS(map_file)
  }
  
  # filter to genes of interest
  annotation.uni <- annotation.uni %>%
    filter(hgnc_symbol %in% genes)
  queryBy <- list(molType = "protein", idType = "UniProt_ID", ids = annotation.uni$uniprotswissprot)
  qresult2 <- drugTargetAnnot(queryBy, config=config)
  qresult2 <- qresult2 %>%
    filter(Action_Type == "INHIBITOR",
           !is.na(First_Approval))
  
  # add gene name back to table
  qresult2 <- qresult2 %>%
    inner_join(annotation.uni, by = c("UniProt_ID" = "uniprotswissprot"))
  
  # merge with differentially expressed genes (up genes in our case)
  qresult2 <- qresult2 %>%
    inner_join(diffexpr_genes, by = c("hgnc_symbol" = "genes")) %>% 
    rename("Comparison" = "comparison",
           "Gene" = "hgnc_symbol",
           "Drug" = "Drug_Name") %>%
    mutate(Source = "FDA") %>%
    filter(logFC > 0) %>%
    unique()
  return(qresult2)
}
