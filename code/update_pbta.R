# Author: Komal S. Rathi
# Date: 02/28/2020
# Function: Code to update pbta histologies file

suppressPackageStartupMessages(library(dotenv))

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
pbta_dir <- file.path(ref_dir, 'pbta')
d3b_toolkit_repo <- file.path(root_dir, '../d3b-analysis-toolkit/scripts')

# load env variables
load_dot_env(file.path(d3b_toolkit_repo, '.envrc'))
host <- Sys.getenv('DWH_HOSTNAME')
usr <- Sys.getenv('DWH_USERNAME')
pwd <- Sys.getenv('DWH_PASSWORD')
db <- Sys.getenv('DWH_DATABASE')
schema <- Sys.getenv('DWH_SCHEMATABLE')
py_script <- file.path(d3b_toolkit_repo, 'select-all-pbta-histologies.py')
system(paste("python", py_script, 
             "-o", file.path(pbta_dir, 'pbta-histologies-base-adapt.tsv'),
             "-u", usr,
             "-p", pwd,
             "-n", host,
             "-d", db,
             "-s", schema, sep = " "))