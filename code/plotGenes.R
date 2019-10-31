#############################
# Function to plot Expression
#############################

plotGenes <- function(myRNASeqAnalysisOut = RNASeqAnalysisOut) {
  geneData <- myRNASeqAnalysisOut[[1]][[2]]
  geneData[,"Direction"] <- ifelse(geneData[,"Z_Score"]>0, "Up", "Down")
  geneData[,"Gene"] <- rownames(geneData)
  geneData <- geneData[order(geneData[,"Z_Score"]),]
  geneData[,"Gene"] <- factor(geneData[,"Gene"], levels=geneData[,"Gene"])
  p <- ggplot(geneData, aes(factor(Gene), y=Z_Score, fill=Direction))+geom_bar(stat="identity")+coord_flip()+theme_bw()
  p <- p+xlab("Gene Symbol")+scale_fill_manual(values = c("forest green", "red"))
  return(p)
}