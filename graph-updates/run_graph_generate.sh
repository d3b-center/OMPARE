#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript -e "rmarkdown::render('01-pathway_enrichment_pediatric.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('02-mutational_analysis_maf_cnv_pediatric.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('03-pathway_enrichment_adult.Rmd', clean = TRUE)"
Rscript -e "rmarkdown::render('04-mutational_analysis_maf_cnv_adult.Rmd', clean = TRUE)"
