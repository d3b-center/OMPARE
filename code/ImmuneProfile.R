########################################################
# calculate immune scores only if doesn't exist already
########################################################

calc.immune.scores <- function(fname){
  # create directory under project directory
  dir <- paste0(topDir,'/ImmuneScores')
  cmd <- paste0("mkdir -p ", dir)
  system(cmd)
  
  # resAll <- res
  # print(head(resAll)[1:4])
  rownames(resAll) <- resAll[,1]
  combGenes <- intersect(rownames(resAll), rownames(expData))
  resAll <- cbind(resAll[combGenes,], expData[combGenes,"FPKM"])
  colnames(resAll)[ncol(resAll)] <- "PatSample"
  resAll[,"max"] <- apply(resAll[3:ncol(resAll)], FUN=max, MARGIN=1)
  resAll <- resAll[order(-resAll[,"max"]),]
  resAll <- resAll[!duplicated(resAll[,2]),]
  rownames(resAll) <- resAll[,2]
  resAll <- resAll[-1:-2]
  resAll <- resAll[-ncol(resAll)]
  # raw.scores <- quietly(rawEnrichmentAnalysis(as.matrix(resAll),
  #                                     xCell.data$signatures,
  #                                     xCell.data$genes))
  # raw.scores[,"CellType"] <- rownames(raw.scores)
  # write.table(raw.scores, fname, quote = F, sep = "\t", row.names = F)
  raw.scores <- capture.output(rawEnrichmentAnalysis(as.matrix(resAll),
                                              xCell.data$signatures,
                                              xCell.data$genes, file.name = fname, 
                                              parallel.type = 'FORK'), file = '/dev/null')
  return(raw.scores)
}

ImmuneProfile <- function() {
  fname <- paste0(topDir,'/ImmuneScores/rawScores.txt')
  # calc.immune.scores(fname)
  # if file does not exist, create one
  if(!file.exists(fname)){
    calc.immune.scores(fname)
  }
  raw.scores <- read.delim(fname, check.names = F)
  raw.scores[,"CellType"] <- raw.scores[,1]
  raw.scores[,1] <- NULL
  raw.scoresTS <- gather(raw.scores, "Sample", "Score", -CellType)
  raw.scoresTS[,"IsSample"] <- ifelse(grepl("PatSample", raw.scoresTS[,"Sample"]), T, F)
  p <- ggplot(raw.scoresTS, aes(CellType, Score))+geom_boxplot(outlier.shape=NA)+theme_bw()+ theme(axis.text.x = element_text(angle = 75, hjust = 1))
  raw.scoresTSSample <- raw.scoresTS[raw.scoresTS[,4]==T,]
  p <- p+geom_point(data=raw.scoresTSSample, aes(CellType, Score), colour = "red", size=3, shape="triangle")
  p <- p+theme(axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14,face="bold"))
  return(p)
}