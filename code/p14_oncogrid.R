# Author: Komal S. Rathi
# Date: 08/10/2020
# Function: Script to generate Oncogrid matrix/additional files using CNV, SNV, Fusion and Expression data
library(readxl)

# read existing oncogrid data
oncogrid.path <- file.path('data', 'Reference', 'oncogrid')
if(!file.exists(file.path(oncogrid.path, 'oncoprint_snv_fusion_mat.txt')) & 
   !file.exists(file.path(oncogrid.path, 'oncoprint_cnv_dge_mat.txt'))){
  n <- FALSE
} else {
  n <- TRUE
}

# read only cohort3 data else if present, then read cohort3 + additional 008 generated data
if(n == FALSE){
  snv_fusion_mat <- read.delim(file.path(oncogrid.path,'cohort3_snv_fusions_mat.txt'), check.names = F)
  cnv_dge_mat <- read.delim(file.path(oncogrid.path,'cohort3_cnv_dge_mat.txt'), check.names = F)
  sym <- ';'
  tmb_mat <- read.delim(file.path(oncogrid.path, 'sample_cohort3.TMB'), check.names = F)
  annot_info <- read.delim(file.path(oncogrid.path, 'sample_cohort3.info'), check.names = F)
} else { 
  snv_fusion_mat <- read.delim(file.path(oncogrid.path,'oncoprint_snv_fusion_mat.txt'), check.names = F)
  cnv_dge_mat <- read.delim(file.path(oncogrid.path,'oncoprint_cnv_dge_mat.txt'), check.names = F)
  sym <- ', '
  tmb_mat <- read.delim(file.path(oncogrid.path, 'tmb_mat.txt'), check.names = F)
  annot_info <- read.delim(file.path(oncogrid.path, 'annot_mat.txt'), check.names = F)
}

# format matrices from wide to long format
snv_fusion_mat <- snv_fusion_mat %>%
  dplyr::rename(Sample = 1) %>%
  mutate_all(list(~na_if(.,""))) %>%
  gather("Gene", "label", -Sample) %>%
  mutate(label = strsplit(as.character(label), sym)) %>% 
  unnest(label) %>%
  unique()

cnv_dge_mat <- cnv_dge_mat %>%
  dplyr::rename(Sample = 1) %>%
  mutate_all(list(~na_if(.,""))) %>%
  gather("Gene", "label", -Sample) %>%
  mutate(label = strsplit(as.character(label), sym)) %>% 
  unnest(label) %>%
  unique()

# read reference gene lists
snv <- read.delim('data/Reference/oncogrid/snv-genes', header = F)
fusion <- read.delim('data/Reference/oncogrid/fusion_genes', header = F)
cnv <- read.delim('data/Reference/oncogrid/copy_number_gene', header = F)
deg <- read.delim('data/Reference/oncogrid/dge_gene_list', header = F)

# fill in details for PNOC008 patient of interest
# 1. get degene info PNOC008 patient vs GTEx Brain
fname <- file.path(topDir, 'Summary', paste0(sampleInfo$subjectID, '_summary.xlsx'))
genes.df.up <- readxl::read_xlsx(path = fname, sheet = "DE_Genes_Up")
genes.df.up <- genes.df.up %>%
  mutate(label = "OVE", 
         Gene = Gene_name) %>%
  filter(Comparison == "GTExBrain_1152",
         Gene_name %in% deg$V1) %>%
  dplyr::select(Gene, label)
genes.df.down <- readxl::read_xlsx(path = fname, sheet = "DE_Genes_Down")
genes.df.down <- genes.df.down %>%
  mutate(label = "UNE",
         Gene = Gene_name) %>%
  filter(Comparison == "GTExBrain_1152",
         Gene_name %in% deg$V1) %>%
  dplyr::select(Gene, label)
deg.genes <- rbind(genes.df.up, genes.df.down)

# 2. get cnv info from cnvGenes
if(nrow(cnvGenes) > 1){
  cnv.genes.gain <- cnvGenes %>%
    filter(Gene %in% cnv$V1, 
           Status == "gain") %>%
    mutate(label = "GAI") %>%
    dplyr::select(Gene, label)
  cnv.genes.loss <- cnvGenes %>%
    filter(Gene %in% cnv$V1) %>%
    mutate(label = "LOS") %>%
    dplyr::select(Gene, label)
  cnv.genes <- rbind(cnv.genes.gain, cnv.genes.loss)
}

# 3. get snv info
if(nrow(mutData) > 0){
  mut.genes <- mutData %>%
    mutate(label = case_when(Variant_Classification %in% "Missense_Mutation" ~ "MIS",
                             Variant_Classification %in% "Nonsense_Mutation" ~ "NOS",
                             Variant_Classification %in% "Nonstop_Mutation" ~ "NOT",
                             Variant_Classification %in% "Frame_Shift_Del" ~ "FSD",
                             Variant_Classification %in% "Frame_Shift_Ins" ~ "FSI",
                             Variant_Classification %in% "In_Frame_Del" ~ "IFD",
                             Variant_Classification %in% "Splice_Site" ~ "SPS"),
           Gene = Hugo_Symbol) %>%
    filter(Hugo_Symbol %in% snv$V1,
           !Variant_Classification %in% c("3'Flank", "5'Flank", "3'UTR", "5'UTR", "IGR", "Intron", "RNA")) %>%
    dplyr::select(Gene, label)
}


# 4. get fusion info
if(nrow(fusData) > 0){
  fus.genes <- data.frame(Gene = unique(c(fusData$HeadGene, fusData$TailGene)))
  fus.genes <- fus.genes %>%
    mutate(label = "FUS") %>%
    filter(Gene %in% fusion$V1)
}

# combine fus + snv
snv_fus <- rbind(mut.genes, fus.genes)
snv_fus <- snv_fus %>%
  mutate(Sample = sampleInfo$subjectID)

# combine deg + cnv
cnv_deg <- rbind(cnv.genes, deg.genes)
cnv_deg <- cnv_deg %>%
  mutate(Sample = sampleInfo$subjectID)
  
# add to existing data
snv_fusion_mat <- unique(rbind(snv_fusion_mat, snv_fus))
cnv_dge_mat <- unique(rbind(cnv_dge_mat, cnv_deg))

# uniquify rows
snv_fusion_mat <- snv_fusion_mat %>%
  group_by(Sample, Gene) %>%
  summarise(label = toString(label))
cnv_dge_mat <- cnv_dge_mat %>%
  group_by(Sample, Gene) %>%
  summarise(label = toString(label))

# convert to matrix
snv_fusion_mat <- snv_fusion_mat %>%
  spread(key = Gene, value = 'label') %>%
  column_to_rownames('Sample')
write.table(x = snv_fusion_mat %>% 
              rownames_to_column('Sample'), file = file.path(oncogrid.path, "oncoprint_snv_fusion_mat.txt"), quote = F, sep = "\t", row.names = TRUE)

cnv_dge_mat <- cnv_dge_mat %>%
  spread(key = Gene, value = 'label') %>%
  column_to_rownames('Sample')
write.table(x = cnv_dge_mat %>% 
              rownames_to_column('Sample'), file = file.path(oncogrid.path, "oncoprint_cnv_dge_mat.txt"), quote = F, sep = "\t", row.names = F)

# now add an * to common genes 
colnames(cnv_dge_mat) <- ifelse(colnames(cnv_dge_mat) %in% colnames(snv_fusion_mat), paste0(colnames(cnv_dge_mat),'*'), colnames(cnv_dge_mat))

# I. now merge both matrices
oncogrid.mat <- snv_fusion_mat %>%
  rownames_to_column('Sample') %>%
  full_join(cnv_dge_mat %>%
              rownames_to_column('Sample'), by = "Sample")
write.table(oncogrid.mat, file = file.path(oncogrid.path, 'oncoprint_mat.txt'), quote = F, sep = "\t", row.names = F)


# II. add TMB info
tmb_mat <- tmb_mat %>% 
  dplyr::rename(Sample = 1)
tmb_mat_p <- data.frame(Sample = sampleInfo$subjectID, TMB = tmb.calculate()/tmb)
tmb_mat <- unique(rbind(tmb_mat, tmb_mat_p))
tmb_mat <- tmb_mat[match(oncogrid.mat$Sample, tmb_mat$Sample),]
write.table(tmb_mat, file = file.path(oncogrid.path, 'tmb_mat.txt'), quote = F, sep = "\t", row.names = F)

# III. add annotation info
annot_info <- annot_info %>% 
  dplyr::rename(Sample = 1)
annot_info_p <- data.frame(Sample = sampleInfo$subjectID,
                           Sequencing_Experiment = "WXS,RNA-Seq",
                           Cohort = "PNOC008",
                           Tumor_Descriptor = "Primary",
                           Integrated_Diagnosis = "High_grade_glioma",
                           OS_Status = "LIVING")
annot_info <- unique(rbind(annot_info, annot_info_p))
annot_info <- annot_info[match(oncogrid.mat$Sample, annot_info$Sample),]
write.table(annot_info, file = file.path(oncogrid.path, 'annot_mat.txt'), quote = F, sep = "\t", row.names = F)
