##############################################################
# Purpose: Main code to generate all tables/figures for report
# Author: Pichai Raman
# Date: 3/21/2019
##############################################################

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

#############################
# Source all code 
#############################
# sourced in order of requirement
source(file.path(utils_dir, 'pubTheme.R'))                      # univeral plot theme
source(file.path(utils_dir, 'helper.R'))                        # universal helper functions like z-score etc
source(file.path(utils_dir, 'createCopyNumber.R'))              # create a map of gene symbol and copy number
source(file.path(utils_dir, 'filterCNV.R'))                     # filter CNV
source(file.path(utils_dir, 'filterFusions.R'))                 # filter Fusions
source(file.path(utils_dir, 'filterMutations.R'))               # filter Mutations
source(file.path(utils_dir, "runRNASeqAnalysis.R"))             # run RNA-seq analysis
source(file.path(utils_dir, 'highConfidenceFindingsTable.R'))   # high confidence findings

# Page 1
source(file.path(code_dir, 'p1_patientSampleInfo.R'))           # patient/sample information
source(file.path(code_dir, 'p1_keyClinicalFindingsTable.R'))    # key findings table
source(file.path(code_dir, 'p1_diseaseSpecificInformation.R'))  # disease specific information
source(file.path(code_dir, 'p1_filterGermline.R'))              # filter and output germline data
source(file.path(code_dir, 'p1_genomicSummary.R'))              # genomic summary table

# Page 2
source(file.path(code_dir, 'p2_plotGenes.R'))                   # plot barchart of top 20 up/down genes
source(file.path(code_dir, 'p2_plotPathway.R'))                 # plot barchat of top 10 up/down pathways (sorted by p-value)

# Page 3
source(file.path(code_dir, 'p3_tumorSignaturePlot.R'))          # plot tumor mutational signatures
if(!is.null(tmb)){                                              # plot tumor mutational burden
  source(file.path(code_dir, 'p3_tmbProfile.R'))
}

# Prepare OpenPBTA and TCGA data to be used for functions in Page 4-7
source(file.path(utils_dir, 'batch_correct.R'))                 # batch correction function
source(file.path(utils_dir, 'pbta_format.R'))                   # format PBTA data
source(file.path(utils_dir, 'tcga_format.R'))                   # format TCGA data

# Page 4 (OpenPBTA) and Page 7 (Genomically Similar PNOC008 and OpenPBTA)
source(file.path(utils_dir, 'ImmuneProfile.R'))                 # plot immune profile using xCell
source(file.path(code_dir, 'p4_TISProfile.R'))                  # plot TIS profile

# Page 5 (PNOC008 + OpenPBTA) and Page 6 (PNOC008 + TCGA GBM)
source(file.path(utils_dir, 'getDimRedPlot.R'))                 # plot UMAP     
source(file.path(utils_dir, 'getSimilarPatients.R'))            # get top 20 correlated patients
source(file.path(utils_dir, 'getKMPlot.R'))                     # plot KM plot of top 20 most correlated patients with survival

# Page 8 (PNOC008)
source(file.path(code_dir, 'p8_tabulate_pathways.R'))           # enriched pathways of PNOC008 patients from top 20 genomically similar

# Page 9 (PNOC008 + OpenPBTA)
source(file.path(code_dir, 'p9_ssGSEA.R'))                      # ssGSEA on top 20 genomically similar patients

# Page 10 
source(file.path(code_dir, 'p10_allFindingsTable.R'))           # all findings table

# Page 11
source(file.path(code_dir, 'p11_plotCNV.R'))
source(file.path(code_dir, 'p11_plotNetwork.R'))
source(file.path(code_dir, 'p11_plotCircos.R'))                 # source and run in driver
if(exists('fusData')){
  plotCircos(topDir = topDir)
}

# Page 12 (source and run in driver)
source(file.path(code_dir, 'p12_cnv_exp_heatmap.R'))
if(!file.exists(file.path(topDir, 'complexHeatmap_phgg.png'))){
  create.heatmap(fname = file.path(topDir, 'complexHeatmap_phgg.png'), 
                 genelist = genelist.heatmap$pHGG_Gene_List, 
                 plot.layout = "h")
}
if(!file.exists(file.path(topDir, 'complexHeatmap_cgs.png'))){
  create.heatmap(fname = file.path(topDir, 'complexHeatmap_cgs.png'), 
                 genelist = genelist.heatmap$Cancer_Gene_Census_CNVs, 
                 plot.layout = "h")
}

# Page 13 Mutational Analysis
source(file.path(code_dir, 'p13_mutational_analysis.R'))

# Page 14 Oncogrid 
if(!file.exists(file.path(topDir, 'complexHeatmap_oncogrid.png'))){
  source(file.path(code_dir, 'p14_oncogrid.R'))
  source(file.path(utils_dir, 'plot_oncogrid.R'))
}
