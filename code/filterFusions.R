####################################################
# Function to filter fusions-
####################################################
filterFusions_star <- function(myFusData=fusData, myCancerGenes=cancerGenes, myJunctionReads=2) {
  fusDataFilt <- myFusData
  
  #FunctionToSplitGenes
  splitMyGene <- function(x)
  {
    x <- strsplit(x, split="\\^")[[1]][1]
    return(x)
  }
  
  #Filter by Number of Reads
  fusDataFilt <- fusDataFilt[fusDataFilt[,"JunctionReads"]>myJunctionReads,]
  
  #Filter by Cancer Gene Census
  myCancerGenes <- as.character(myCancerGenes[,1])
  fusDataFilt[,"HeadGene"] <- sapply(as.character(fusDataFilt[,"LeftGene"]), FUN=splitMyGene)
  fusDataFilt[,"TailGene"] <- sapply(as.character(fusDataFilt[,"RightGene"]), FUN=splitMyGene)
  fusDataFilt <- fusDataFilt[fusDataFilt[,"HeadGene"]%in%myCancerGenes|fusDataFilt[,"TailGene"]%in%myCancerGenes,]
  
  return(fusDataFilt)
  
}

filterFusions_aribba <- function(myFusData=fusData, myCancerGenes=cancerGenes)
{
  fusDataFilt <- myFusData
  fusDataFilt[,"X.fusion_name"] <- paste(as.character(fusDataFilt[,"X.gene1"]), "-", as.character(fusDataFilt[,"gene2"]), sep="")
  fusDataFilt[,"Splice_type"] <- fusDataFilt[,"type"]
  
  
  #Filter by Number of Reads
  fusDataFilt <- fusDataFilt[fusDataFilt[,"confidence"]!="low",]
  
  #Filter by Cancer Gene Census
  myCancerGenes <- as.character(myCancerGenes[,1])
  fusDataFilt[,"HeadGene"] <- as.character(fusDataFilt[,"X.gene1"])
  fusDataFilt[,"TailGene"] <- as.character(fusDataFilt[,"gene2"])
  fusDataFilt <- fusDataFilt[fusDataFilt[,"HeadGene"]%in%myCancerGenes|fusDataFilt[,"TailGene"]%in%myCancerGenes,]
  
  return(fusDataFilt)
  
}

####################################################
# End Function to filter fusions
####################################################