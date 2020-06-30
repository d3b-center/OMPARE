##############################################################
# Purpose: Main code to generate all tables/figures for report
# Author: Pichai Raman
# Date: 3/21/2019
##############################################################

#############################
# Source all code 
#############################
# sourced in order of requirement
source('code/pubTheme.R')                       # univeral plot theme
source('code/helper.R')                         # universal helper functions like z-score etc
source('code/patientSampleInfo.R')              # functions to load patient info
source('code/createCopyNumber.R')               # create a map of gene symbol and copy number
source('code/filterDruggability.R')             # create a map of gene symbol and drug info from dgiDB
source('code/filterCNV.R')                      # filter CNV
source('code/filterFusions.R')                  # filter Fusions
source('code/filterMutations.R')                # filter Mutations
source("code/runRNASeqAnalysis.R")              # function to run RNA-seq analysis
RNASeqAnalysisOut <- runRNASeqAnalysis(expData) # run in driver
source('code/highConfidenceFindingsTable.R')    # high confidence findings

# Page 1
source('code/p1_patientSampleInfo.R')           # patient/sample information
source('code/p1_keyClinicalFindingsTable.R')    # key findings table
source('code/p1_diseaseSpecificInformation.R')  # disease specific information
source('code/p1_filterGermline.R')              # filter and output germline data
source('code/p1_genomicSummary.R')              # genomic summary table

# Page 2
source('code/p2_plotGenes.R')                   # plot barchart of top 20 up/down genes
source('code/p2_plotPathway.R')                 # plot barchat of top 10 up/down pathways (sorted by p-value)

# Page 3
source('code/p3_tumorSignaturePlot.R')          # plot tumor mutational signatures
if(!is.null(tmb)){                              # plot tumor mutational burden
  source('code/p3_tmbProfile.R')
}

# Prepare OpenPBTA and TCGA data to be used for functions in Page 4-7
source('code/pbta_format.R')                    # format PBTA data
source('code/tcga_format.R')                    # format TCGA data

# Page 4 (OpenPBTA) and Page 7 (Genomically Similar PNOC008 and OpenPBTA)
source('code/ImmuneProfile.R')                  # plot immune profile using xCell
source('code/p4_TISProfile.R')                  # plot TIS profile

# Page 5 (PNOC008 + OpenPBTA) and Page 6 (PNOC008 + TCGA GBM)
source('code/getTSNEPlot.R')                    # plot t-SNE     
source('code/getSimilarPatients.R')             # get top 20 correlated patients
source('code/getKMPlot.R')                      # plot KM plot of top 20 most correlated patients with survival

# Page 8 (PNOC008)
source('code/p8_tabulate_pathways.R')           # enriched pathways of PNOC008 patients from top 20 genomically similar

# Page 9 (PNOC008 + OpenPBTA)
source('code/p9_ssGSEA.R')                      # ssGSEA on top 20 genomically similar patients

# Page 10 
source('code/p10_allFindingsTable.R')           # all findings table

# Page 11
source('code/p11_plotCNV.R')
source('code/p11_plotNetwork.R')
source('code/p11_plotCircos.R')                 # source and run in driver
if(exists('fusData')){
  plotCircos(topDir = topDir)
}

# Page 12 (source and run in driver)
source('code/p12_cnv_exp_heatmap.R')
if(!file.exists(paste0(topDir,'/complexHeatmap_phgg.png'))){
  create.heatmap(fname = paste0(topDir,'/complexHeatmap_phgg.png'), 
                 genelist = genelist.heatmap$pHGG_Gene_List, 
                 plot.layout = "h")
}
if(!file.exists(paste0(topDir,'/complexHeatmap_cgs.png'))){
  create.heatmap(fname = paste0(topDir,'/complexHeatmap_cgs.png'), 
                 genelist = genelist.heatmap$Cancer_Gene_Census_CNVs, 
                 plot.layout = "h")
}

# Page 13 Mutational Analysis
source('code/p13_mutational_analysis.R')