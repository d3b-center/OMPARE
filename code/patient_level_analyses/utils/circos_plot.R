circos_plot <- function(topDir = topDir, chrMap, cancerGenes, fname) {
  
  # filter chrMap to chr 1:22, X and Y
  chrMap.subset <- chrMap %>% 
    filter(chromosome %in% c(1:22,'X','Y'),
           hgnc_symbol != "") %>%
    mutate(chromosome = paste0("chr", chromosome))
  
  # set up core components  
  data(UCSC.HG38.Human.CytoBandIdeogram) 
  RCircos.Set.Core.Components(cyto.info = UCSC.HG38.Human.CytoBandIdeogram, 
                              chr.exclude = NULL, 
                              tracks.inside = 10, 
                              tracks.outside = 0)
  
  png(filename = fname, height = 8, width = 9, units = "in", res = 300)
  # initialize plot
  RCircos.Set.Plot.Area()  
  RCircos.Chromosome.Ideogram.Plot()
  
  # 1. add mutations
  mut.genes <- unique(mutDataFilt$Hugo_Symbol)
  if(length(mut.genes) > 0){
    genomic.data <- chrMap.subset %>% 
      filter(hgnc_symbol %in% mut.genes) %>%
      dplyr::select(chromosome, gene_start, gene_end, hgnc_symbol) %>%
      arrange(hgnc_symbol)
    RCircos.Gene.Connector.Plot(genomic.data = genomic.data, track.num = 1, side = "in")
    RCircos.Gene.Name.Plot(gene.data = genomic.data, name.col = 4, track.num = 2, side = "in")
  }
  
  # 2. add expression 
  expr.data <- data.frame('hgnc_symbol' = names(rnaseq_analysis_output$expr.genes.logfc), 
                          'expression' = rnaseq_analysis_output$expr.genes.logfc)
  expr.data$expression <- ifelse(expr.data$expression > 5, 5, ifelse(expr.data$expression < (-5), -5, expr.data$expression))
  
  # create heatmap with all genes
  heatmap.data <- chrMap.subset %>%
    inner_join(expr.data, by = "hgnc_symbol") %>%
    dplyr::select(chromosome, gene_start, gene_end, hgnc_symbol, expression)
  RCircos.Heatmap.Plot(heatmap.data = heatmap.data, data.col = 5, track.num = 4, side = "in")
  
  # only label highly differential (abs(z-score) > 5)
  genomic.data <- heatmap.data %>%
    filter(abs(expression) == 5,
           hgnc_symbol %in% cancerGenes$Gene_Symbol)
  RCircos.Gene.Connector.Plot(genomic.data = genomic.data, track.num = 5, side = "in")
  RCircos.Gene.Name.Plot(gene.data = genomic.data, name.col = 4, track.num = 6, side = "in")
  
  # 3. add fusions 
  if(exists('fusData')){
    fusion.data <- fusData %>%
      filter(!grepl(",", TailGene),
             !grepl(",", HeadGene)) %>%
      filter(HeadGene %in% chrMap.subset$hgnc_symbol | TailGene %in% chrMap.subset$hgnc_symbol) %>%
      unique()
    
    if(nrow(fusion.data) > 0){
      fusion.data <- fusion.data %>%
        inner_join(chrMap.subset, by = c("HeadGene" = "hgnc_symbol")) %>%
        inner_join(chrMap.subset, by = c("TailGene" = "hgnc_symbol")) %>%
        dplyr::select(X.fusion_name, 
                      HeadGene, gene_start.x, gene_end.x, chromosome.x,
                      TailGene, gene_start.y, gene_end.y, chromosome.y) 
      
      # add link between fusion coordinates
      link.data <- fusion.data %>%
        dplyr::select(chromosome.x, gene_start.x, gene_end.x,
                      chromosome.y, gene_start.y, gene_end.y)
      RCircos.Link.Plot(link.data = link.data, track.num = 12, by.chromosome = TRUE)
      
      # add fusion gene names
      gene.data.head <- fusion.data %>%
        dplyr::select(chromosome.x, gene_start.x, gene_end.x, HeadGene) %>%
        setNames(., c("chromosome", "gene_start", "gene_end", "hgnc_symbol"))
      gene.data.tail <- fusion.data %>%
        dplyr::select(chromosome.y, gene_start.y, gene_end.y, TailGene) %>%
        setNames(., c("chromosome", "gene_start", "gene_end", "hgnc_symbol"))
      gene.data <- rbind(gene.data.head, gene.data.tail)
      RCircos.Gene.Name.Plot(gene.data = gene.data, name.col = 4, track.num = 9, inside.pos = 50)
    }
  }
  
  dev.off()
}