########################################################
# calculate immune scores only if doesn't exist already
########################################################

calc_immune_scores <- function(fullmat, fname){
  # create directory under project directory
  raw.scores <- capture.output(rawEnrichmentAnalysis(as.matrix(fullmat),
                                              xCell.data$signatures,
                                              xCell.data$genes, 
                                              file.name = fname, 
                                              parallel.type = 'FORK'), file = '/dev/null')
  return(raw.scores)
}

immune_profile <- function(fullmat, fname) {
  # if file does not exist, create one
  if(!file.exists(fname)){
    calc_immune_scores(fullmat, fname)
  }
  raw.scores <- read.delim(fname, check.names = F)
  raw.scores <- raw.scores %>%
    dplyr::rename(CellType = 1) %>%
    gather("Sample", "Score", -CellType) %>%
    mutate("IsSample" = ifelse(grepl(sampleInfo$subjectID, Sample), T, F))
  raw.scores.sample <- raw.scores %>%
    filter(IsSample == TRUE)
  
  # set factors
  celltype.order <- raw.scores %>%
    group_by(CellType) %>%
    summarise(median = median(Score)) %>%
    arrange(desc(median)) %>%
    .$CellType
  raw.scores$CellType <- factor(raw.scores$CellType, levels = celltype.order)
  
  # boxplot
  p <- ggplot(raw.scores, aes(CellType, Score)) + 
    geom_boxplot(outlier.shape = NA) +  
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    geom_point(data = raw.scores.sample, aes(CellType, Score), colour = "red", size = 3, shape = "triangle") +
    theme(axis.text = element_text(size = 8, face = "bold"), 
          axis.title = element_blank())
  return(p)
}
