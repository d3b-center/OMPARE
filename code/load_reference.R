#########################
# Load all reference data
#########################

# GTEx (7859 samples)
# Read GTEx Expression Data
gtexData <- readRDS("data/Reference/GTEx/GTEx_fullExpr_matrix.RDS")
# Read GTEx Annotation
gtexGeneAnnot <- read.delim("data/Reference/GTEx/gencode.v23.annotation.gi_ti_gs.txt", stringsAsFactors =F)
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
load("data/Reference/cbttc_genes_fpkm_1110.RData")
res[,1] <- sapply(as.character(res[,1]), FUN=remDotStuff)
mapping <- read.delim("data/Reference/mappingFile.txt", header=F, stringsAsFactors=F)
clinData <- read.delim("data/Reference/study_view_clinical_data.txt", stringsAsFactors=F)
tmbData <- read.csv("data/reference/complete_results.csv", stringsAsFactors=F)
if(!file.exists('data/Reference/survData.txt')){
  formatData <- parseSurvival()
  write.table(formatData, "data/Reference/survData.txt", row.names=F, sep="\t")
} 
survData <- read.delim("data/Reference/survData.txt", stringsAsFactors=F)
signatures <- readAlexandrovSignatures("data/Reference/signatures_probabilities.txt")
dgidb <- read.delim("data/Reference/DGIdb.txt", stringsAsFactors = F)
rawSurvData <- read.delim("data/Reference/CBTTC_PullApril21.txt", stringsAsFactors = F)
# germlineMarkers <- read.delim("data/Reference/germlineMarkers.txt", stringsAsFactors = F)
pharmacogenomics.genes <- read.delim("data/Reference/Pharmacogenomics_Genes.list", stringsAsFactors = F, header = F)
pharmacogenomics.genes <- data.frame("Gene" = pharmacogenomics.genes$V1, "Class" = "Pharmacogenomics")
chop.panel.genes <- read.delim("data/Reference/CHOP_Additional_Cancer_Genes.list", stringsAsFactors = F, header = F)
chop.panel.genes <-  data.frame("Gene" = chop.panel.genes$V1, "Class" = "CHOP Panel")
acmg.genes <- read.delim("data/Reference/ACMG_22_Cancer_Genes.list", stringsAsFactors = F, header = F)
acmg.genes <-  data.frame("Gene" = acmg.genes$V1, "Class" = "ACMG")
germlineMarkers <- rbind(pharmacogenomics.genes,  acmg.genes, chop.panel.genes)
