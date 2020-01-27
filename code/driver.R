##############################################################
# Purpose: Main code to generate all tables/figures for report
# Author: Pichai Raman
# Date: 3/21/2019
##############################################################


#############################
# Load packages 
#############################
suppressPackageStartupMessages(library(flexdashboard))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(copynumber))
suppressPackageStartupMessages(library(RCircos))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(GeneNetworkBuilder))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(network))
suppressPackageStartupMessages(library(sna))
suppressPackageStartupMessages(library(ggnetwork))
suppressPackageStartupMessages(library(Rtsne))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(randomcoloR))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(RTCGA.clinical))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(xCell))
suppressPackageStartupMessages(library(decompTumor2Sig))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(GenomicRanges))

#############################
# Source all code 
#############################
# sourced in order of requirement
# source('code/helper.R')                         # load misc functions
source('code/createCopyNumber.R')               
source('code/filterDruggability.R')
source('code/parseSurvival.R')
source('code/germlineAnalysis.R')
source('code/patientSampleInfo.R')              # load patient info
source('code/pubTheme.R')                       # plot theme
# source('code/load_reference.R')                 # load all reference data
source('code/filterCNV.R')                      # filter CNV
source('code/filterFusions.R')                  # filter Fusions
source('code/filterMutations.R')                # filter Mutations
source("code/runRNASeqAnalysis.R")              # RNA-seq analysis


#############################
# Code Blocks
#############################
# Function to get RNA-Seq and Pathway Analysis
if(exists('expData')){
  RNASeqAnalysisOut <- runRNASeqAnalysis(expData) # *Run in driver
}
source('code/plotGenes.R')
source('code/plotPathway.R')

# Key Clinical Findings P1
source('code/key_findings.R')

# High Confidence Alterations P2
source('code/highConfidenceFindingsTable.R')

# Immmune Signatures P3
source('code/ImmuneProfile.R')

# TMB Tumor Signatures P4
source('code/tmbProfile.R')
source('code/tumorSignaturePlot.R')

# Genomically Similar Samples P5
if(exists('expData')){
  source('code/cleanData_p5.R')                   # *Run in driver
}
source('code/getTSNEPlot.R')
source('code/getKMPlot.R')
source('code/getSimilarPatients.R')

# All findings table P6
source('code/allFindingsTable.R')

# Genomic Landscape P7
source('code/genomic_landscape.R')

# only run if fusion data is present
if(exists('fusData')){
  plotCircos(topDir = topDir)
}

# only run if fusion or mutation data is present
if(exists('expData')){
  if(exists('mutData') | exists('fusData')){
    plotNetwork()
  }
}
