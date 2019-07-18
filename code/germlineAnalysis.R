############################
# Purpose: Germline Analysis
############################

germlineAnalysis <- function(mutData, germlineMarkers){
  # germlineData <- read.delim("~/Projects/OMPARE/data/PNOC008/PNOC008-2/a9ecf1e7-a7ac-441b-ac37-a2930966600f.gatk_vqsr.vep.maf", skip=1, stringsAsFactors = F);
  # germlineData <- read.delim("~/Projects/OMPARE/data/CBTTC-HGG/MutationsMAF/5c3eace5-950a-4a05-81ed-5c04b4a0a367.strelka.vep.maf", skip = 1, stringsAsFactors = F)
  germlineData <- mutData
  
  # Filter to only genes in germline markers
  germlineMarkersGenes <- unique(germlineMarkers[,"Gene"]);
  germlineData <- germlineData[germlineData[,1] %in% germlineMarkersGenes,]
  
  # Filter by annotation
  germlineData <- germlineData[germlineData[,"Consequence"]!="synonymous_variant",]
  germlineData <- germlineData[grepl("pathogenic", germlineData[,"CLIN_SIG"]),]
  germlineOut <- merge(germlineMarkers, germlineData, by.y ="Hugo_Symbol", by.x="Gene")
  
  germlineOut <- germlineOut[,c("Disease", "Gene", "HGVSp_Short")]
  return(germlineOut)
}
