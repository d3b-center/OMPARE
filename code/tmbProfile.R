tmbProfile <- function(myTMB=tmbData) {
  
  #Determine if sample is an outlier in dataset
  isOutlier <- function(distTmp) {
    c25 <- as.numeric(quantile(distTmp, .25))
    c75 <- as.numeric(quantile(distTmp, .75))
    myThresh <- (c75-c25)+c75
    out <- distTmp>myThresh
  }
  
  
  getMMRFreq <- function(tmpAnnot) {
    #Define the Hypermutation profile
    tmpAnnot[,"Outlier"] <- isOutlier(tmpAnnot[,"MUTATION_COUNT"])
    tmpAnnot[,"EST_MUT_PER_MB"] <- tmpAnnot[,"MUTATION_COUNT"]/30
    return(tmpAnnot)
  }
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  
  ######################
  ##Read data 
  ####################
  
  #Read in data 26504 to start
  
  ######################
  ##Filter data 
  ####################
  
  #Filter to studies that have mutation counts 15373
  myTMB[,"MUTATION_COUNT"] <- as.numeric(myTMB[,"MUTATION_COUNT"])
  myTMB <- myTMB[!is.na(myTMB[,"MUTATION_COUNT"]),]
  myTMB[,"NAME"] <- as.character(myTMB[,"NAME"])
  
  #Remove PPTP and FMI
  myTMB <- filter(myTMB, DISEASE.NAME!="Pediatric Preclinical Testing Program" & DISEASE.NAME!="Mixed Cancer Types")
  
  #Add Whether Pediatric or Adult
  myTMB[,"Type"] <- "Adult"
  pediatricCancers <- c("Neuroblastoma", "Medulloblastoma", "Acute Myeloid Leukemia", "Pediatric Ganglio", "Diffuse Intrinsic Pontine Glioma")
  myTMB[myTMB[,"DISEASE.NAME"]%in%pediatricCancers, "Type"] <- "Pediatric"
  
  #Sort and order it in decreasing fashion 
  grouped <- group_by(myTMB, Type, DISEASE.NAME)
  tmpOut <- summarise(grouped, median=median(MUTATION_COUNT))
  tmpOut <- arrange(tmpOut, Type, desc(median))
  myOrder <- as.character(as.data.frame(tmpOut)[,"DISEASE.NAME"])
  myTMB[,"DISEASE.NAME"] <- factor(myTMB[,"DISEASE.NAME"], levels=myOrder)
  myTMB[,"EST_MUT_PER_MB"] <- myTMB[,"MUTATION_COUNT"]/30
  
  #Filter mutations to get TMB of samples
  filtMut <- mutData 
  filtMut <- filtMut[filtMut[,"HGVSp_Short"]!="",]
  tmbSample <- nrow(filtMut)/30
  
  #Plot it
  p <- ggplot(myTMB, aes(DISEASE.NAME, EST_MUT_PER_MB, fill=Type))+geom_boxplot()+theme_bw()
  p <- p+scale_y_log10(breaks=c(.25, 1, 10, 100, 500))+scale_fill_manual(values=c("blue", "red"))
  p <- p+xlab("Disease")+ylab("Est Mutation Count per MB (Log Scale)")+theme(axis.text.x = element_text(angle = -90, hjust = (0)))
  p <- p+ggtitle("Mutation Count versus Disease")+geom_hline(yintercept=tmbSample, linetype=2)
  return(p)
  
}