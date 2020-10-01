########################################################
# calculate immune scores only if doesn't exist already
########################################################

calc.immune.scores <- function(fullmat, fname){
  # create directory under project directory
  raw.scores <- capture.output(rawEnrichmentAnalysis(as.matrix(fullmat),
                                              xCell.data$signatures,
                                              xCell.data$genes, file.name = fname, 
                                              parallel.type = 'FORK'), file = '/dev/null')
  return(raw.scores)
}

ImmuneProfile <- function(fullmat, fname) {
  # fname <- paste0(topDir,'/ImmuneScores/rawScores.txt')
  # calc.immune.scores(fname)
  # if file does not exist, create one
  if(!file.exists(fname)){
    calc.immune.scores(fullmat, fname)
  }
  raw.scores <- read.delim(fname, check.names = F)
  raw.scores[,"CellType"] <- raw.scores[,1]
  raw.scores[,1] <- NULL
  raw.scoresTS <- gather(raw.scores, "Sample", "Score", -CellType)
  raw.scoresTS[,"IsSample"] <- ifelse(grepl(sampleInfo$subjectID, raw.scoresTS[,"Sample"]), T, F)
  p <- ggplot(raw.scoresTS, aes(CellType, Score)) + 
    geom_boxplot(outlier.shape = NA) +  
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 75, hjust = 1))
  raw.scoresTSSample <- raw.scoresTS[raw.scoresTS$IsSample == T,]
  p <- p + 
    geom_point(data = raw.scoresTSSample, aes(CellType, Score), colour = "red", size = 3, shape = "triangle") +
    theme(axis.text = element_text(size = 8, face = "bold"), 
          axis.title = element_blank())
  return(p)
}
