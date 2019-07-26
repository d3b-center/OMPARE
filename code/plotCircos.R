##########################
# Function for Circos Plot
##########################

plotCircos <- function(topDir = topDir) {
  data(UCSC.HG38.Human.CytoBandIdeogram) 
  chr.exclude <- NULL
  cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
  tracks.inside <- 10
  tracks.outside <- 0
  RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
  
  rcircos.params <- RCircos.Get.Plot.Parameters()
  rcircos.cyto <- RCircos.Get.Plot.Ideogram()
  rcircos.position <- RCircos.Get.Plot.Positions()
  
  fname <- paste0(topDir, "/tmpRCircos.png")
  png(fname)
  RCircos.Set.Plot.Area()  
  #  par(mai=c(0.25, 0.25, 0.25, 0.25))
  plot.new()
  plot.window(c(-2.5,2.5), c(-2.5, 2.5))
  RCircos.Chromosome.Ideogram.Plot()
  
  # Put in Mutations
  mySymbols=as.character(filterMutations()[,1])
  RCircos.Gene.Label.Data <- chrMap[chrMap[,1]%in%mySymbols,]
  RCircos.Gene.Label.Data <- RCircos.Gene.Label.Data[!grepl("CHR_", RCircos.Gene.Label.Data[,4]), ]
  RCircos.Gene.Label.Data[,4] <- paste("chr", RCircos.Gene.Label.Data[,4], sep="")
  colnames(RCircos.Gene.Label.Data) <- c("hgnc_symbol", "start_position", "end_position", "chromosome_name")
  RCircos.Gene.Label.Data <- RCircos.Gene.Label.Data[,c("chromosome_name", "start_position", "end_position", "hgnc_symbol")]
  RCircos.Gene.Label.Data <- RCircos.Gene.Label.Data[order(RCircos.Gene.Label.Data[,4]),]
  name.col <- 4
  side <- "in"
  track.num <- 1
  RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
  track.num <- 2
  RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col,track.num, side)
  
  # Add Expression
  tmpExp <- data.frame(names(RNASeqAnalysisOut[[1]][[1]]), RNASeqAnalysisOut[[1]][[1]])
  tmpExp[,2] <- ifelse(tmpExp[,2]>5, 5, ifelse(tmpExp[,2]<(-5), -5, tmpExp[,2]))
  colnames(tmpExp) <- c("GeneName", "Expression")
  RCircos.Heatmap.Data <- merge(chrMap, tmpExp, by.x="HGNC.symbol", by.y="GeneName")
  colnames(RCircos.Heatmap.Data) <- c("GeneName", "chromStart", "chromEnd", "Chromosome", "Expression")
  RCircos.Heatmap.Data <- RCircos.Heatmap.Data[!grepl("CHR_", RCircos.Heatmap.Data[,4]), ]
  RCircos.Heatmap.Data <- RCircos.Heatmap.Data[!grepl("MT", RCircos.Heatmap.Data[,4]), ]
  RCircos.Heatmap.Data[,4] <- paste("chr", RCircos.Heatmap.Data[,4], sep="")
  RCircos.Heatmap.Data <- RCircos.Heatmap.Data[,c("Chromosome", "chromStart", "chromEnd", "GeneName", "Expression")]
  data.col <- 5
  track.num <- 4
  side <- "in"
  RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, track.num, side)
  RCircos.Heatmap.Data.High <- RCircos.Heatmap.Data[abs(RCircos.Heatmap.Data[,"Expression"])==5,]
  RCircos.Heatmap.Data.High <- RCircos.Heatmap.Data.High[RCircos.Heatmap.Data.High[,"GeneName"]%in%as.character(cancerGenes[,1]),]
  RCircos.Gene.Connector.Plot(RCircos.Heatmap.Data.High, 5, side)
  RCircos.Gene.Name.Plot(RCircos.Heatmap.Data.High, 4,6, side)
  
  # Add Fusions
  myFus <- fusData
  # fusData <- fusData[which(fusData$HeadGene %in% chrMap$HGNC.symbol & fusData$TailGene %in% chrMap$HGNC.symbol),]
  RCircos.Link.Data.tmp.h <- chrMap[chrMap[,1]%in%myFus[,"HeadGene"],]
  RCircos.Link.Data.tmp.h <- RCircos.Link.Data.tmp.h[!grepl("CHR_", RCircos.Link.Data.tmp.h[,"Chromosome.scaffold.name"]),]
  RCircos.Link.Data.tmp.h <- RCircos.Link.Data.tmp.h[,c(4,2,3,1)]
  RCircos.Link.Data.tmp.h[,1] <- paste("chr", RCircos.Link.Data.tmp.h[,1], sep="")
  
  RCircos.Link.Data.tmp.t <- chrMap[chrMap[,1]%in%myFus[,"TailGene"],]
  RCircos.Link.Data.tmp.t <- RCircos.Link.Data.tmp.t[!grepl("CHR_", RCircos.Link.Data.tmp.t[,"Chromosome.scaffold.name"]),]
  RCircos.Link.Data.tmp.t <- RCircos.Link.Data.tmp.t[,c(4,2,3,1)]
  RCircos.Link.Data.tmp.t[,1] <- paste("chr", RCircos.Link.Data.tmp.t[,1], sep="")
  
  RCircos.Link.Data <- data.frame(RCircos.Link.Data.tmp.h[,c(1,2,3)], RCircos.Link.Data.tmp.t[,c(1,2,3)])
  track.num <- 12
  RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE)
  RCircos.Gene.Name.Plot(rbind(RCircos.Link.Data.tmp.h, RCircos.Link.Data.tmp.t), 4,9, inside.pos=50)
  dev.off()
}