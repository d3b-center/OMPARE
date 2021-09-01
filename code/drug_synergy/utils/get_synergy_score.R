

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("igraph"))

#### Define Directories ----------------------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

analysis_dir <- file.path(root_dir, "code", "drug_synergy")
results_dir <- file.path(analysis_dir, "results", "synergy_score")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}

#### Functions Used ----------------------------------------------------------

#' Title Drug synergy score calculation function (for calculating synergy score of drug 1 and drug 2)
#'
#' @param tar1 Gene targets of drug 1 that are in the subnetwork
#' @param tar2 Gene targets of drug 2 that are in the subnetwork
#' @param sigNet1 This is the network induced using subnetwork of the module (the module that is positively correlated with the cluster of our sample of interest)
#'                Subnetwork is induced using `graph.edgelist` function from igraph package.
#'
#' @return sctar1 is the sum of centrality scores of all gene targets of drug 1 in subnetwork
#'         sctar1 is the sum of centrality scores of all gene targets of drug 2 in subnetwork
#'         sScore is the sum of sctar1 and sctar2
#' @export 
#'
#' @examples getSynScore2(targets_in_sub1, targets_in_sub2, subnetwork_graphed) 
#'

getSynScore2 <- function(tar1, tar2, sigNet1){
  
  n0 <- length(tar1);
  n1 <- length(tar2);
  vSet0 <- V(sigNet1)$name
  
  ## calculate target centralization score using closeness and betweenness 
  sc1 <- closeness(sigNet1, vSet0);  sc2 <- betweenness(sigNet1, vSet0, directed=F); sc3 <- page.rank(sigNet1, algo="prpack", vids=vSet0)$vector; ### for normalization
  sc11 <- closeness(sigNet1, tar1);  sc12 <- betweenness(sigNet1, tar1, directed=F); sc13 <- page.rank(sigNet1, algo="prpack", vids=tar1)$vector;
  sc21 <- closeness(sigNet1, tar2);  sc22 <- betweenness(sigNet1, tar2, directed=F); sc23 <- page.rank(sigNet1, algo="prpack", vids=tar2)$vector;
  sc11 <- (sc11-min(sc1))/(max(sc1)-min(sc1)); sc12 <- (sc12-min(sc2))/(max(sc2)-min(sc2)); sc13 <- (sc13-min(sc3))/(max(sc3)-min(sc3));          ###normalization
  sc21 <- (sc21-min(sc1))/(max(sc1)-min(sc1)); sc22 <- (sc22-min(sc2))/(max(sc2)-min(sc2)); sc23 <- (sc23-min(sc3))/(max(sc3)-min(sc3));
  sctar1 <- (sum(sc11)+sum(sc12)+sum(sc13))/3.0;
  sctar2 <- (sum(sc21)+sum(sc22)+sum(sc23))/3.0;
  
  sScore <- sctar1 + sctar2;
  return(sScore);
  return(sctar1);
  return(sctar2);
  
}

