suppressPackageStartupMessages({
  library(RCircos)
})

circos_plot <- function(chr_map, cancer_genes, fname, filtered_mutations, rnaseq_analysis_output, filtered_fusions) {
  
  # filter chr_map to chr 1:22, X and Y
  chr_map_subset <- chr_map %>% 
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
  mut_genes <- unique(filtered_mutations$Hugo_Symbol)
  if(length(mut_genes) > 0){
    genomic_data <- chr_map_subset %>% 
      filter(hgnc_symbol %in% mut_genes) %>%
      dplyr::select(chromosome, gene_start, gene_end, hgnc_symbol) %>%
      arrange(hgnc_symbol)
    RCircos.Gene.Connector.Plot(genomic.data = genomic_data, track.num = 1, side = "in")
    RCircos.Gene.Name.Plot(gene.data = genomic_data, name.col = 4, track.num = 2, side = "in")
  }
  
  # 2. add expression 
  expr_data <- data.frame('hgnc_symbol' = names(rnaseq_analysis_output$expr.genes.logfc), 
                          'expression' = rnaseq_analysis_output$expr.genes.logfc)
  expr_data$expression <- ifelse(expr_data$expression > 5, 5, ifelse(expr_data$expression < (-5), -5, expr_data$expression))
  
  # create heatmap with all genes
  heatmap_data <- chr_map_subset %>%
    inner_join(expr_data, by = "hgnc_symbol") %>%
    dplyr::select(chromosome, gene_start, gene_end, hgnc_symbol, expression)
  RCircos.Heatmap.Plot(heatmap.data = heatmap_data, data.col = 5, track.num = 4, side = "in")
  
  # only label highly differential (abs(z-score) > 5)
  genomic_data <- heatmap_data %>%
    filter(abs(expression) == 5,
           hgnc_symbol %in% cancer_genes$Gene_Symbol)
  RCircos.Gene.Connector.Plot(genomic.data = genomic_data, track.num = 5, side = "in")
  RCircos.Gene.Name.Plot(gene.data = genomic_data, name.col = 4, track.num = 6, side = "in")
  
  # 3. add fusions 
  if(exists('filtered_fusions')){
    fusion_data <- filtered_fusions %>%
      filter(!grepl(",", gene2),
             !grepl(",", gene1)) %>%
      filter(gene1 %in% chr_map_subset$hgnc_symbol | gene2 %in% chr_map_subset$hgnc_symbol) %>%
      unique()
    
    if(nrow(fusion_data) > 0){
      fusion_data <- fusion_data %>%
        inner_join(chr_map_subset, by = c("gene1" = "hgnc_symbol")) %>%
        inner_join(chr_map_subset, by = c("gene2" = "hgnc_symbol")) %>%
        dplyr::select(fusion_name, 
                      gene1, gene_start.x, gene_end.x, chromosome.x,
                      gene2, gene_start.y, gene_end.y, chromosome.y) 
      
      # add link between fusion coordinates
      link_data <- fusion_data %>%
        dplyr::select(chromosome.x, gene_start.x, gene_end.x,
                      chromosome.y, gene_start.y, gene_end.y) %>% 
        as.data.frame()
      
      if(nrow(link_data) > 0){
        RCircos.Link.Plot(link.data = link_data, track.num = 12, by.chromosome = TRUE)
        
        # add fusion gene names
        gene_data_head <- fusion_data %>%
          dplyr::select(chromosome.x, gene_start.x, gene_end.x, gene1) %>%
          setNames(., c("chromosome", "gene_start", "gene_end", "hgnc_symbol"))
        gene_data_tail <- fusion_data %>%
          dplyr::select(chromosome.y, gene_start.y, gene_end.y, gene2) %>%
          setNames(., c("chromosome", "gene_start", "gene_end", "hgnc_symbol"))
        gene_data <- rbind(gene_data_head, gene_data_tail) %>% as.data.frame()
        RCircos.Gene.Name.Plot(gene.data = gene_data, name.col = 4, track.num = 9, inside.pos = 50)
      }
    }
  }
  
  dev.off()
}
