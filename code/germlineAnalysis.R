###########################################
#Purpose: Germline Analysis
#Author: Pichai Raman
#Date: 3/21/2019
###########################################


germlineData <- read.delim("../data/PNOC008/PNOC008-2/a9ecf1e7-a7ac-441b-ac37-a2930966600f.gatk_vqsr.vep.maf", skip=1, stringsAsFactors = F);

germlineMarkers <- read.delim("../data/Reference/germlineMarkers.txt", stringsAsFactors = F);


#Filter to only genes in germline markers
germlineMarkersGenes <- unique(germlineMarkers[,"Gene"]);
germlineData <- germlineData[germlineData[,1]%in%germlineMarkersGenes,]

#Filter by annotation
germlineData <- germlineData[germlineData[,"Consequence"]!="synonymous_variant",]
germlineData <- germlineData[grepl("pathogenic", germlineData[,"CLIN_SIG"]),]
germlineOut <- merge(germlineMarkers, germlineData, by.y ="Hugo_Symbol", by.x="Gene")


germlineOut <- germlineOut[,c("Disease", "Gene", "HGVSp_Short")]