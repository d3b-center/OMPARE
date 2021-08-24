# PBTA adapt histology file
# remove unwanted samples
pbta_hist_adapt <- read.delim('../data/pbta-histologies-base-adapt.tsv')
pbta_hist_adapt <- pbta_hist_adapt %>%
  filter(kf_visibility == "True") %>%
  dplyr::select(-c(tumor_fraction, tumor_ploidy))

pbta_hist <- read.delim('../data/pbta-histologies.tsv')
pbta_hist <- pbta_hist %>%
  dplyr::select(Kids_First_Biospecimen_ID, tumor_fraction) %>%
  unique() %>%
  left_join(pbta_hist_adapt, by = 'Kids_First_Biospecimen_ID')

# subset to cancer genes 
cancer_genes <- readRDS('../data/cancer_gene_list.rds')
gene.list <- unique(cancer_genes$Gene_Symbol)

# read controlfreec (add tumor purity from histology file and get tumor ploidy from controlfreec)
pbta.controlfreec <- data.table::fread('../data/pbta-cnv-controlfreec.tsv.gz')
pbta.controlfreec <- pbta.controlfreec %>%
  dplyr::select(Kids_First_Biospecimen_ID, tumor_ploidy) %>%
  unique() %>%
  inner_join(pbta_hist, by = "Kids_First_Biospecimen_ID")

# read cnv kit and combine with above
pbta.cnvkit <- data.table::fread('../data/pbta-cnv-cnvkit.seg.gz')
pbta.cnvkit <- pbta.cnvkit %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "ID") %>%
  inner_join(pbta.controlfreec, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::rename("SampleID" = "sample_id")

# chr coordinates to gene symbol map
chr_map <- read.delim('../data/mart_export_genechr_mapping.txt', stringsAsFactors = F, check.names = F)
colnames(chr_map) <- c("hgnc_symbol", "gene_start", "gene_end", "chromosome")
chr_map <- chr_map %>%
  filter(hgnc_symbol != "") %>%
  mutate(chromosome = paste0("chr", chromosome))

# overlap
library(GenomicRanges)
subject <- with(chr_map, GRanges(chromosome, IRanges(start = gene_start, end = gene_end, names = hgnc_symbol)))
query <- with(pbta.cnvkit, GRanges(chrom, IRanges(start = loc.start, end = loc.end)))

# find overlaps and subset maf 
res <- findOverlaps(query = query, subject = subject, type = "within")
pbta.cnvkit <- data.frame(pbta.cnvkit[queryHits(res),], chr_map[subjectHits(res),])

filter_cnv <- function(myCNVData, myCancerGenes = cancer_genes) {
  
  # tumor suppressor genes
  myTSGenes <- myCancerGenes %>%
    filter(grepl(pattern = "TumorSuppressorGene|Is.Tumor.Suppressor.Gene", type)) %>%
    .$Gene_Symbol %>%
    unique()
  
  # oncogenes
  myOncogenes <- myCancerGenes %>%
    filter(grepl(pattern = "Oncogene|OncoKB.Annotated|Is.Oncogene", type)) %>%
    .$Gene_Symbol %>%
    unique()
  myOncogenes <- setdiff(myOncogenes, myTSGenes)
  
  # gain in oncogenes
  cnvDataFiltUp <- myCNVData %>%
    filter(status %in% c("Gain", "Amplification") & hgnc_symbol %in% myOncogenes)
  
  # loss in tsgs
  cnvDataFiltDown <- myCNVData %>%
    filter(status %in% c("Loss", "Complete Loss") & hgnc_symbol %in% myTSGenes)
  
  # combine
  cnvDataFilt <- rbind(cnvDataFiltUp, cnvDataFiltDown)
  
  return(cnvDataFilt)
}

# cnvData <- pbta.cnvkit %>% filter(Kids_First_Biospecimen_ID == "BS_01Y5F4PN")
merge_cnv <- function(cnvData, cancer_genes){
  sample_name <- unique(cnvData$Kids_First_Biospecimen_ID)
  print(sample_name)
  ploidy <- cnvData %>% .$tumor_ploidy %>% unique() %>% as.numeric()
  
  # add copy number status
  cnvData <- cnvData %>%
    mutate(status = case_when(copy.num == 0 ~ "Complete Loss",
                              copy.num < ploidy & copy.num > 0 ~ "Loss",
                              copy.num == ploidy ~ "Neutral",
                              copy.num > ploidy & copy.num < ploidy + 3 ~ "Gain",
                              copy.num >= ploidy + 3 ~ "Amplification"))
  
  cnvDataFilt <- filter_cnv(myCNVData = cnvData, myCancerGenes = cancer_genes)
  
  # add sample name
  if(nrow(cnvDataFilt) > 0) {
    cnvDataFilt$sample_name <- sample_name
  }
  
  return(cnvDataFilt)
}

pbta.cnv <- plyr::ddply(.data = pbta.cnvkit, 
                        .variables = 'Kids_First_Biospecimen_ID', 
                        .fun = function(x) merge_cnv(cnvData = x, cancer_genes = cancer_genes))
pbta.cnv <- pbta.cnv %>%
  mutate(Alteration_Datatype = "CNV",
         Alteration_Type = stringr::str_to_title(status),
         Alteration = paste0('Copy Number Value:', copy.num),
         Study = "PBTA",
         Gene = hgnc_symbol) %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  unique()
saveRDS(pbta.cnv, file = '../data/pbta-cnv-cnvkit-filtered.rds')
