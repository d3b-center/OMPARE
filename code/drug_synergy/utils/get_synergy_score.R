

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

