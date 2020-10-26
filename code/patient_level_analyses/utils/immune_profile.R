########################################################
# calculate immune scores only if doesn't exist already
########################################################

source(file.path(patient_level_analyses_utils, 'quiet.R'))

immune_profile <- function(fullmat, fname) {
  # immune scores
  if(!file.exists(fname)){
    raw.scores <- as.data.frame(quiet(xCellAnalysis(expr = fullmat)))
    write.table(raw.scores, file = fname, quote = F, sep = "\t")
  } else {
    raw.scores <- read.delim(fname, check.names = F)
  }
  
  # format data
  raw.scores <- raw.scores %>%
    rownames_to_column('CellType') %>%
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
