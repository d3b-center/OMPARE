# this script is a modified version of drugTargetInteractions::getParalogs

getParalogs_modified <- function (queryBy) {
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  IDMresult <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                             "uniprot_gn_symbol", "uniprotswissprot", "uniprotsptrembl", 
                                             "description"), filters = queryBy$idType, values = queryBy$ids, mart = mart)
  IDMresult <- data.frame(QueryID = as.character(IDMresult$ensembl_gene_id), 
                          ENSEMBL = as.character(IDMresult$ensembl_gene_id), GENES = as.character(IDMresult$uniprot_gn_symbol), 
                          ID_up_sp = as.character(IDMresult$uniprotswissprot), 
                          ID_up_sp_tr = as.character(IDMresult$uniprotsptrembl))
  result <- biomaRt::getBM(attributes = c("external_gene_name", 
                                          "ensembl_gene_id", "hsapiens_paralog_associated_gene_name", 
                                          "hsapiens_paralog_ensembl_gene", "hsapiens_paralog_perc_id", 
                                          "hsapiens_paralog_perc_id_r1"), filters = queryBy$idType, 
                           values = queryBy$ids, mart = mart)
  query <- result[!duplicated(result$ensembl_gene_id), ]
  query <- data.frame(QueryID = as.character(query$ensembl_gene_id), 
                      ENSEMBL = as.character(query$ensembl_gene_id), GENES = as.character(query$external_gene_name), 
                      stringsAsFactors = FALSE)
  query <- cbind(query, paralog_perc_id = 100, paralog_perc_id_r1 = 100)
  
  # result$hsapiens_paralog_ensembl_gene[is.na(result$hsapiens_paralog_ensembl_gene)] <- "" # added by KSR
  result <- result[nchar(result$hsapiens_paralog_ensembl_gene) >  0, ]
  paralog <- data.frame(QueryID = result$ensembl_gene_id, ENSEMBL = result$hsapiens_paralog_ensembl_gene, 
                        GENES = result$hsapiens_paralog_associated_gene_name, 
                        paralog_perc_id = result$hsapiens_paralog_perc_id, paralog_perc_id_r1 = result$hsapiens_paralog_perc_id_r1, 
                        stringsAsFactors = FALSE)
  resultDF <- rbind(query, paralog)
  resultDF <- resultDF[order(resultDF$QueryID), ]
  tmp <- resultDF[complete.cases(resultDF), ]
  uniprot <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                           "uniprot_gn_symbol", "uniprotswissprot", "uniprotsptrembl", 
                                           "description"), filters = "ensembl_gene_id", 
                            values = resultDF$ENSEMBL, 
                            mart = mart)
  # added by KSR
  resultDF <- resultDF %>%
    as.data.frame() %>%
    dplyr::filter(ENSEMBL %in% uniprot$ensembl_gene_id)
  up_sp <- tapply(uniprot$uniprotswissprot, uniprot$ensembl_gene_id, 
                  function(x) unique(as.character(x)[nchar(x) > 0]), simplify = FALSE)
  index <- vapply(up_sp, length, integer(1))
  up_sp[index == 0] <- ""
  index[index == 0] <- 1
  index <- index[as.character(resultDF$ENSEMBL)]
  up_sp <- up_sp[names(index)]
  resultDF <- cbind(resultDF[rep(seq_along(resultDF$ENSEMBL), 
                                 index), ], ID_up_sp = as.character(unlist(up_sp)))
  up_sp_tr <- tapply(uniprot$uniprotsptrembl, uniprot$ensembl_gene_id, 
                     function(x) unique(as.character(x)[nchar(x) > 0]), simplify = FALSE)
  index <- vapply(up_sp_tr, length, integer(1))
  up_sp_tr[index == 0] <- ""
  index[index == 0] <- 1
  index <- index[as.character(resultDF$ENSEMBL)]
  up_sp_tr <- up_sp_tr[names(index)]
  resultDF <- cbind(resultDF[rep(seq_along(resultDF$ENSEMBL), 
                                 index), ], ID_up_sp_tr = as.character(unlist(up_sp_tr)))
  row.names(resultDF) <- NULL
  IDMresult[IDMresult == ""] <- NA
  for (i in colnames(IDMresult)) {
    if (is.factor(IDMresult[, i])) {
      IDMresult[, i] <- as.character(IDMresult[, i])
    }
  }
  resultDF[resultDF == ""] <- NA
  for (i in colnames(resultDF)) {
    if (is.factor(resultDF[, i])) {
      resultDF[, i] <- as.character(resultDF[, i])
    }
  }
  resultList <- list(IDM = IDMresult, SSNN = resultDF)
  return(resultList)
}
