# Prepare PBTA data for Mutational analysis
# This needs to be run with every updated version of OpenPBTA (current v21)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "filter_cnv.R"))
pbta.dir <- file.path(root_dir, 'data', 'pbta')

# PBTA adapt histology file
# remove unwanted samples
pbta_hist_adapt <- read.delim(file.path(pbta.dir, 'pbta-histologies-base-adapt.tsv'))
pbta_hist_adapt <- pbta_hist_adapt %>%
  filter(kf_visibility == "True") %>%
  dplyr::select(-c(tumor_fraction, tumor_ploidy))

pbta_hist <- read.delim(file.path(pbta.dir, 'pbta-histologies.tsv'))
pbta_hist <- pbta_hist %>%
  dplyr::select(Kids_First_Biospecimen_ID, tumor_fraction) %>%
  unique() %>%
  left_join(pbta_hist_adapt, by = 'Kids_First_Biospecimen_ID')

# subset to cancer genes 
cancer_genes <- readRDS(file.path(root_dir, "data", "cancer_gene_list.rds"))
gene.list <- unique(cancer_genes$Gene_Symbol)

# read controlfreec (add tumor purity from histology file and get tumor ploidy from controlfreec)
pbta.controlfreec <- data.table::fread(file.path(pbta.dir, 'pbta-cnv-controlfreec.tsv.gz'))
pbta.controlfreec <- pbta.controlfreec %>%
  dplyr::select(Kids_First_Biospecimen_ID, tumor_ploidy) %>%
  unique() %>%
  inner_join(pbta_hist, by = "Kids_First_Biospecimen_ID")

# read cnv kit and combine with above
pbta.cnvkit <- data.table::fread(file.path(pbta.dir, 'pbta-cnv-cnvkit.seg.gz'))
pbta.cnvkit <- pbta.cnvkit %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "ID") %>%
  inner_join(pbta.controlfreec, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::rename("SampleID" = "sample_id")

# chr coordinates to gene symbol map
chr_map <- read.delim(file.path(root_dir, "data", 'mart_export_genechr_mapping.txt'), stringsAsFactors = F, check.names = F)
colnames(chr_map) <- c("hgnc_symbol", "gene_start", "gene_end", "chromosome")
chr_map <- chr_map %>%
  filter(hgnc_symbol != "") %>%
  mutate(chromosome = paste0("chr", chromosome))

# overlap
subject <- with(chr_map, GRanges(chromosome, IRanges(start = gene_start, end = gene_end, names = hgnc_symbol)))
query <- with(pbta.cnvkit, GRanges(chrom, IRanges(start = loc.start, end = loc.end)))

# find overlaps and subset maf 
res <- findOverlaps(query = query, subject = subject, type = "within")
pbta.cnvkit <- data.frame(pbta.cnvkit[queryHits(res),], chr_map[subjectHits(res),])

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
saveRDS(pbta.cnv, file = file.path(pbta.dir, 'pbta-cnv-cnvkit-filtered.rds'))

# read fusions (oncogenic arriba fusions)
pbta.fusions <- read.delim(file.path(pbta.dir, 'pbta-fusion-putative-oncogenic.tsv'))
pbta.fusions <- pbta.fusions[grep('ARRIBA', pbta.fusions$CalledBy),]
pbta.fusions <- pbta.fusions  %>%
  unite(col = "Gene", Gene1A, Gene1B, Gene2A, Gene2B, sep = ", ", na.rm = T)  %>%
  mutate(Alteration_Datatype = "Fusion",
         Alteration_Type = stringr::str_to_title(Fusion_Type),
         Alteration = FusionName,
         Kids_First_Biospecimen_ID = Sample,
         Study = "PBTA") %>%
  inner_join(pbta_hist, by = "Kids_First_Biospecimen_ID") %>%
  mutate(SampleID = sample_id) %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  separate_rows(Gene, convert = TRUE) %>%
  filter(Gene %in% gene.list) %>%
  unique()
saveRDS(pbta.fusions, file = file.path(pbta.dir, 'pbta-fusion-putative-oncogenic-filtered.rds'))

# read mutations (consensus)
pbta.mutations <- fread(file.path(pbta.dir, 'pbta-snv-consensus-mutation.maf.tsv.gz'))
keepVC <- c("Nonsense_Mutation", "Missense_Mutation", 
            "Splice_Region", "Splice_Site",
            "3'UTR", "5'UTR", 
            "In_Frame_Ins", "In_Frame_Del",
            "Frame_Shift_Ins", "Frame_Shift_Del")
keepVI <- c("MODIFIER", "MODERATE", "HIGH")
pbta.mutations <- pbta.mutations %>%
  filter(BIOTYPE == "protein_coding",
         Variant_Classification %in% keepVC,
         IMPACT %in% keepVI,
         Hugo_Symbol %in% gene.list) %>%
  mutate(Gene = Hugo_Symbol,
         Alteration_Datatype = "Mutation",
         Alteration_Type = Variant_Classification,
         Alteration = HGVSp_Short,
         Kids_First_Biospecimen_ID = Tumor_Sample_Barcode,
         Study = "PBTA") %>%
  inner_join(pbta_hist, by = "Kids_First_Biospecimen_ID") %>%
  mutate(SampleID = sample_id) %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  unique()
saveRDS(pbta.mutations, file = file.path(pbta.dir, 'pbta-snv-consensus-mutation-filtered.rds'))

