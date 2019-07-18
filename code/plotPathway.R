##########################
# Function to plot pathway
##########################

plotPathway <- function(myRNASeqAnalysisOut=RNASeqAnalysisOut) {
  pathData <- myRNASeqAnalysisOut[[2]][[2]]
  pathDataUp <- pathData[pathData[,"Direction"]=="Up",]
  pathDataUp[,"Pathway"] <- rownames(pathDataUp)
  pathDataUp <- pathDataUp[order(pathDataUp[,"P_VAL"]),]
  pathDataUp[,"Pathway"] <- factor(pathDataUp[,"Pathway"], levels=pathDataUp[,"Pathway"])
  pathDataUp <- pathDataUp[1:10,]
  
  pathDataDown <- pathData[pathData[,"Direction"]=="Down",]
  pathDataDown[,"Pathway"] <- rownames(pathDataDown)
  pathDataDown <- pathDataDown[order(pathDataDown[,"P_VAL"]),]
  pathDataDown[,"Pathway"] <- factor(pathDataDown[,"Pathway"], levels=pathDataDown[,"Pathway"])
  pathDataDown <- pathDataDown[1:10,]
  pathData <- rbind(pathDataDown, pathDataUp)
  pathData[,"Direction"] <- factor(pathData[,"Direction"], levels=c("Down", "Up"))
  p <- ggplot(pathData, aes(factor(Pathway), y=(-1)*log10(P_VAL), fill=Direction))+geom_bar(stat="identity")+coord_flip()+theme_bw()
  p <- p+xlab("Pathway Name")+scale_fill_manual(values = c("forest green", "red"))+ylab("-log10 P-Value")
  return(p)
}