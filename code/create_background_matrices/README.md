#### Background matrices

Author: [Komal Rathi](https://github.com/komalsrathi)
Code: [create_background_matrices](../create_background_matrices)

##### Introduction 

The master genomics files are saved under `../../data/master_genomics` which contains all pediatric tumor data. 	

```
../../data/master_genomics
├── consensus_wgs_plus_cnvkit_wxs.BS_GQFPB8F3.tsv.gz
├── gene-counts-rsem-expected_count.BS_S3MET0JZ.rds
├── gene-expression-rsem-tpm.BS_S3MET0JZ.rds
├── pbta-BS_S3MET0JZ-fusion-arriba.tsv.gz
├── pbta-BS_S3MET0JZ-fusion-starfusion.tsv.gz
├── pbta_0928_histologies-ops-BS_GQFPB8F3.tsv
└── snv-consensus-plus-hotspots.BS_GQFPB8F3.maf.gz
```

Depending on the input `patient_cancer_type`, a mapping file (`../../data/tumor_normal_mapping.tsv`) is read and the reference `pediatric tumor` type, `adult tumor` type and `normal tissue` is determined. These values are then used inside of `../create_background_matrices` to create the correponding subset of comparator matrices (i.e. expression/mutations/cnv/fusions) for adult tumors, pediatric tumors and normal tissues from `../../data/master_genomics` files and saved under `../../data/adult_data`, `../../data/pediatric_data` and `../../data/normal_data`, respectively. 

The patient of interest data would be included within the pediatric tumor subset and can be extracted by reading the files under `../../data/pediatric_data/cancer_type`. This is done because the input genomic files are growing and reading the full set is time consuming every time a module is called.  

Example for `LGG` and `HGG`:

```
tree ../../data/
adult_data
├── GBM
│   ├── GBM_cnv_filtered.rds
│   ├── GBM_counts.rds
│   ├── GBM_histologies.tsv
│   ├── GBM_mutation_filtered.rds
│   └── GBM_tpm.rds
└── LGG
    ├── LGG_cnv_filtered.rds
    ├── LGG_counts.rds
    ├── LGG_histologies.tsv
    ├── LGG_mutation_filtered.rds
    └── LGG_tpm.rds
normal_data
└── Brain
    ├── Brain_counts.rds
    ├── Brain_histologies.tsv
    └── Brain_tpm.rds
pediatric_data
├── HGG
│   ├── HGG_cnv_filtered.rds
│   ├── HGG_counts.rds
│   ├── HGG_degs.rds
│   ├── HGG_fusion_filtered.rds
│   ├── HGG_histologies.tsv
│   ├── HGG_mutation_filtered.rds
│   └── HGG_tpm.rds
└── LGG
    ├── LGG_cnv_filtered.rds
    ├── LGG_counts.rds
    ├── LGG_degs.rds
    ├── LGG_fusion_filtered.rds
    ├── LGG_histologies.tsv
    ├── LGG_mutation_filtered.rds
    └── LGG_tpm.rds
```

Also, we are utilizing the RNA-seq and histology information for `GTEx` and `TCGA` cohorts from `OpenPedCan-analysis`, so it is being soft linked under `../../data/OpenPedCan-analysis`.

```
../../data/OpenPedCan-analysis
├── analyses -> /Users/rathik/Projects/PediatricOpenTargets/OpenPedCan-analysis/analyses
└── data -> /Users/rathik/Projects/PediatricOpenTargets/OpenPedCan-analysis/data
```

##### Output details

* Clinical (`histologies.tsv`): cancer-type specific histologies 
* Expression (`tpm.rds`/`counts.rds`): cancer-type specific expression matrices collapsed to unique gene symbols and subsetted to protein coding genes.
* Copy number (`cnv_filtered.rds`): cancer-type specific consensus cnv files filtered using [filter_cnv.R](../utils/filter_cnv.R) script.
* Mutation (`mutation_filtered.rds`): cancer-type specific mutation files filtered using [filter_mutations.R](../utils/filter_mutations.R) script.
* Fusions (`fusion_filtered.rds`): cancer-type specific fusion files filtered using [filter_fusions.R](../utils/filter_fusions.R) script.
