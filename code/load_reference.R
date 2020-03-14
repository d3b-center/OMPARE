#########################
# Load all reference data
#########################

# GTEx (7859 samples)
gtexData <- readRDS("data/Reference/GTEx/GTEx_matrix.RDS")

# Cancer Genes
cancerGenes <- read.delim("data/Reference/CancerGeneList.tsv", stringsAsFactors = F)
cancerGenes <- cancerGenes %>%
  filter(Gene_Symbol != "") %>%
  dplyr::select(-Count) %>%
  gather(key = "file", value = "type", -Gene_Symbol) %>%
  mutate(type = file)
geneListRef <- read.delim("data/Reference/genelistreference.txt", stringsAsFactors = F)
geneListRef <- subset(geneListRef, type == "TumorSuppressorGene" | type == "CosmicCensus" | type == "Oncogene")
cancerGenes <- rbind(cancerGenes, geneListRef)
rm(geneListRef)
tsgGenes <- read.delim("data/Reference/Human_TSGs.txt", stringsAsFactors = F)

chrMap <- read.delim("data/Reference/mart_export_genechr_mapping.txt", stringsAsFactors =F)
geneMania <- read.delim("data/Reference/GeneManiaNetwork.txt", stringsAsFactors =F)
diseaseSpecificFields <- read.delim("data/Reference/DiseaseSpecificFields.txt")
hallMarkSets <- getGmt("data/Reference/mSigDB/h.all.v6.2.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
hallMarkSets <- geneIds(hallMarkSets)
hallMarkSetsTS <- stack(hallMarkSets)
res <- readRDS('data/Reference/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds')
clinData <- read.delim("data/Reference/pbta-histologies.tsv", stringsAsFactors = F)
survData <- clinData %>% 
  filter(experimental_strategy == "RNA-Seq") %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, OS_days, OS_status) %>%
  filter(!is.na(OS_status)) %>%
  mutate(OS_status = ifelse(OS_status == "DECEASED", 1, 0))
signatures <- readAlexandrovSignatures("data/Reference/signatures_probabilities.txt")
dgidb <- read.delim("data/Reference/DGIdb.txt", stringsAsFactors = F)

# Germline
pharmacogenomics.genes <- read.delim("data/Reference/Pharmacogenomics_Genes.list", stringsAsFactors = F, header = F)
pharmacogenomics.genes <- data.frame("Gene" = pharmacogenomics.genes$V1, "Class" = "Pharmacogenomics")
chop.panel.genes <- read.delim("data/Reference/CHOP_Additional_Cancer_Genes.list", stringsAsFactors = F, header = F)
chop.panel.genes <-  data.frame("Gene" = chop.panel.genes$V1, "Class" = "CHOP Panel")
acmg.genes <- read.delim("data/Reference/ACMG_22_Cancer_Genes.list", stringsAsFactors = F, header = F)
acmg.genes <-  data.frame("Gene" = acmg.genes$V1, "Class" = "ACMG")
germlineMarkers <- rbind(pharmacogenomics.genes,  acmg.genes, chop.panel.genes)

# TMB
pedTMB <- data.table::fread("data/Reference/pbta-TMBscores_withdiseastype.txt")
adultTMB <- data.table::fread("data/Reference/TCGA_diseasetypes_and_samples_TMBscores.txt")
TMBFileBED <- data.table::fread("data/Reference/xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed")
colnames(TMBFileBED)  <- c("chr", "start", "end")
