## Plot Conversions For Tables in Mutational Analysis Module of OMPARE

Currently, `Mutational Analysis` under both `Pediatric Tumors Analysis` and `Adult Tumors Analysis` are using tabular display. We want to leverage the bioconductor package `maftools` (https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html) to display these results as graphs.

### Usage
After all the relevant files are downloaded and stored in the `data` folder of the main repository, the module can be run with:

```
bash run_graph-generate.sh
```

All data files required are:
- For the patient of interest (download from Cavatica delivery project):
  - `mutational_analysis_adult.rds`
  - `mutational_analysis_pediatric.rds`
  - `pathway_analysis_adult.rds`
  - `pathway_analysis_pediatric.rds`

- Combined file for PNOC008 patients (download from s3://d3b-bix-dev-data-bucket/PNOC008/reference/pnoc008):
  - `pnoc008_clinical.rds`
  - `pnoc008_cnv_filtered.rds`
  - `pnoc008_consensus_maf_combined.rds`

- Combined file for PNOC008 patients (download from s3://d3b-bix-dev-data-bucket/PNOC008/reference/pbta):
  - `pbta-histologies.tsv`
  - `pbta-cnv-cnvkit-filtered.rds`
  - `pbta-snv-consensus-mutation.maf.tsv.gz`

- Combined files for TCGA samples (download from s3://d3b-bix-dev-data-bucket/PNOC008/reference/tcga):
  - `tcga_gbm_clinical.rds`
  - `tcga_gbm_cnv_filtered.rds`

###### Contents

`01-pathway_enrichment_pediatric.Rmd` reads in `pathway_analysis_pediatric.rds` and output `pathway_analysis_pediatric.png` in `plots/pediatric/pathway_barplot` - this is the pathway barplot with only sample of interest. 

`02-mutational_analysis_maf_cnv_pediatric.Rmd` use the samples listed in `mutational_analysis_pediatric.rds`.
`pbta-histologies.tsv` and `pnoc008_clinical.rds` were used for sample id matching.
To get complete MAF results, `pnoc008_consensus_maf_combined.rds` and `pbta-snv-consensus-mutation.maf.tsv.gz` were queried.
To get complete CNV results, `pnoc008_cnv_filtered.rds` and `pbta-cnv-cnvkit-filtered.rds` were queried. 
**NOTE**: `mutational_analysis_pediatric.rds` results are a combination of 20 transcriptomically similar samples as compared to our sample of interest.

For results in both `recurrent_alterations` and `shared_genes` with `Alteration_Type == "Mutation"`, the following plots are generated:
  - Summary plots for MAF:
    - `mutational_recurrent_pediatric.pdf`
    - `mutational_shared_pediatric.pdf`
  - Lollipop plots for top 10 most mutated genes:
    - `lollipop_recurrent_pediatric.pdf`
    - `lollipop_shared_pediatric.pdf`

For results in both `recurrent_alterations` and `shared_genes` with `Alteration_Type == "CNV"`, the following plots are generated:
  - Oncoplots for MAF+CNV:
    - `mutational_cnv_recurrent_pediatric.pdf`
    - `mutational_cnv_shared_pediatric.pdf`

Input:
  - `../data/mutational_analysis_pediatric.rds`
  - `../data/pnoc008_clinical.rds`
  - `../data/pnoc008_cnv_filtered.rds`
  - `../data/pnoc008_consensus_maf_combined.rds`
  - `../data/pbta-histologies.tsv`
  - `../data/pbta-cnv-cnvkit-filtered.rds`
  - `../data/pbta-snv-consensus-mutation.maf.tsv.gz`
Output:
  - `plots/pediatric/lollipop/lollipop_recurrent_pediatric.pdf`
  - `plots/pediatric/lollipop/lollipop_shared_pediatric.pdf`
  - `plots/pediatric/oncoplot/mutational_cnv_recurrent_pediatric.pdf`
  - `plots/pediatric/oncoplot/mutational_cnv_shared_pediatric.pdf`
  - `plots/pediatric/plotmaf/mutational_recurrent_pediatric.pdf`
  - `plots/pediatric/plotmaf/mutational_shared_pediatric.pdf`

`03-pathway_enrichment_adult.Rmd` reads in `pathway_analysis_adult.rds` and output `pathway_analysis_adult.png` in `plots/adult/pathway_barplot` - this is the pathway barplot with only sample of interest. 

`04-mutational_analysis_maf_cnv_adult.Rmd` use the samples listed in `mutational_analysis_adult.rds`.
`tcga_gbm_clinical.rds` and `pnoc008_clinical.rds` were used for sample id matching.
To get complete MAF results, `pnoc008_consensus_maf_combined.rds` were queried. For TCGA samples, mutect2 MAF files were downloaded and modified with each run.
To get complete CNV results, `pnoc008_cnv_filtered.rds` and `tcga_gbm_cnv_filtered.rds` were queried. 
**NOTE**: `mutational_analysis_adult.rds` results are a combination of 20 transcriptomically similar samples as compared to our sample of interest.

For results in both `recurrent_alterations` and `shared_genes` with `Alteration_Type == "Mutation"`, the following plots are generated:
  - Summary plots for MAF:
    - `mutational_recurrent_adult.pdf`
    - `mutational_shared_adult.pdf`
  - Lollipop plots for top 10 most mutated genes:
    - `lollipop_recurrent_adult.pdf`
    - `lollipop_shared_adult.pdf`

For results in both `recurrent_alterations` and `shared_genes` with `Alteration_Type == "CNV"`, the following plots are generated:
  - Oncoplots for MAF+CNV:
    - `mutational_cnv_recurrent_adult.pdf`
    - `mutational_cnv_shared_adult.pdf`

Input:
  - `../data/mutational_analysis_adult.rds`
  - `../data/pnoc008_clinical.rds`
  - `../data/pnoc008_cnv_filtered.rds`
  - `../data/pnoc008_consensus_maf_combined.rds`
  - `../data/tcga_gbm_clinical.rds`
  - `../data/tcga_gbm_cnv_filtered.rds`
Output:
  - `plots/adult/lollipop/lollipop_recurrent_adult.pdf`
  - `plots/adult/lollipop/lollipop_shared_adult.pdf`
  - `plots/adult/oncoplot/mutational_cnv_recurrent_adult.pdf`
  - `plots/adult/oncoplot/mutational_cnv_shared_adult.pdf`
  - `plots/adult/plotmaf/mutational_recurrent_adult.pdf`
  - `plots/adult/plotmaf/mutational_shared_adult.pdf`

`pnoc008_maf_combine.R` is used to generated `pnoc008_consensus_maf_combined.rds` and maintain the file - it is not part of the graph generation pipeline and can be incorporated in current OMPARE. 

`cnvkit_modify.R` is used to generated `pbta-cnv-cnvkit-filtered.rds` and it is adapted from current OMPARE (already incorporated) - just kept here for reference.

