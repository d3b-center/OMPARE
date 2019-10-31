#########################
# Load all reference data
#########################

# Read GTEx Expression Data
gtexData <- readRDS("data/Reference/GTEx/GTEx_fullExpr_matrix.RDS")
# Read GTEx Annotation
gtexGeneAnnot <- read.delim("data/Reference/GTEx/gencode.v23.annotation.gi_ti_gs.txt", stringsAsFactors =F)
cancerGenes <- read.delim("data/Reference/CancerGeneList.tsv", stringsAsFactors = F)
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
germlineMarkers <- read.delim("data/Reference/germlineMarkers.txt", stringsAsFactors = F)
