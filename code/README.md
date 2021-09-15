## Pre-requisites 

1. [Install R packages](https://github.com/d3b-center/OMPARE/blob/master/code/utils/install_pkgs.R)

2. Download files from PNOC008 data delivery project on cavatica with each new patient and save under `OMPARE/results/PNOC008-xx`:

```
    results/PNOC008-xx
    ├── copy-number-variations
    │   ├──.controlfreec.CNVs.p.value.txt
    │   ├──.controlfreec.info.txt
    │   ├──.controlfreec.ratio.txt
    │   ├──.diagram.pdf
    │   └──.gainloss.txt
    ├── gene-expressions
    │   └──.rsem.genes.results.gz (TPM + Expected counts)
    ├── gene-fusions
    │   ├──.STAR.fusion_predictions.abridged.coding_effect.tsv
    │   └──.arriba.fusions.tsv
    └── simple-variants
        ├──.lancet_somatic.norm.annot.protected.maf
        ├──.mutect2_somatic.norm.annot.protected.maf
        ├──.strelka2_somatic.norm.annot.protected.maf
        ├──.vardict_somatic.norm.annot.protected.maf
        ├──.consensus_somatic.protected.maf
        └──.gatk.PASS.vcf.gz.hg38_multianno.txt.gz
```

3. Download PNOC008 Clinical Manifest from the Kids First data tracker to get `PNOC008 subject ID` (i.e. PNOC008-XX), `Research ID` (i.e. cohort_participant_id), `Age at Collection in days` (i.e. OS days) and `Last Known Status` (i.e. OS_status). This file has free text in some of the fields including `diagnosis` so currently I am hardcoding the `short_histology` as `HGAT` and `broad_histology` as `Diffuse astrocytic and oligodendroglial tumor`. 
4. Download PNOC008 Sample Manifest from the Kids First data tracker to get the `library name` (to code for stranded/poly-a needed for batch correction when combining PBTA + PNOC008 or TCGA + PNOC008 matrices for transcriptomically similar and downstream analyses).
5. Download `pbta-histologies-base-adapt.tsv` to map cohort_participant_id from sample/clinical manifests to `Kids_First_Biospecimen_ID`
6. [Merge and update combined matrices for PNOC008 patients](https://github.com/d3b-center/OMPARE/blob/master/code/update_pnoc008_matrices.R). The matrices are updated under `OMPARE/data/pnoc008`:
```
    data/pnoc008
    ├── pnoc008_clinical.rds
    ├── pnoc008_cnv_filtered.rds
    ├── pnoc008_consensus_mutation.rds
    ├── pnoc008_consensus_mutation_filtered.rds
    ├── pnoc008_counts_matrix.rds
    ├── pnoc008_fpkm_matrix.rds
    ├── pnoc008_fusions_filtered.rds
    ├── pnoc008_tmb_scores.rds
    ├── pnoc008_tpm_matrix.rds
    └── pnoc008_vs_gtex_brain_degs.rds
```

## Analysis modules + files needed

[Driver script](https://github.com/d3b-center/OMPARE/blob/master/code/driver.R) is used to load all libraries and reference files as well as call the analyses modules in a specific order.

|                                                                       Analysis                                                                      | Expression  | Fusions | CNV | Mutations |
|:---------------------------------------------------------------------------------------------------------------------------------------------------:|:-----------:|:-------:|:---:|:---------:|
| [Merge + Update PNOC008 matrices with each new patient](https://github.com/d3b-center/OMPARE/blob/master/code/update_pnoc008_matrices.R)                                                                                               |      Y      |    Y    |  Y  |     Y     |
|  [RNAseq analysis: Differential expression + GSEA enrichment + Barplots](https://github.com/d3b-center/OMPARE/tree/master/code/rnaseq_analysis)                                                                              |      Y      |    N    |  N  |     N     |
|  [All/Key clinical findings + Germline variants](https://github.com/d3b-center/OMPARE/tree/master/code/p1_modules)                                                                                                      |      Y      |    Y    |  Y  |     Y     |
|  [Tumor mutational signatures and mutational burden](https://github.com/d3b-center/OMPARE/tree/master/code/tmb_analysis)                                                                                                  |      N      |    N    |  N  |     Y     |
|  [Transcriptomically similar analysis](https://github.com/d3b-center/OMPARE/tree/master/code/transcriptomically_similar_analysis) (UMAP, ssGSEA)                                                                                                  |      Y      |    N    |  N  |     N     |
|  [Transcriptomically similar analysis](https://github.com/d3b-center/OMPARE/tree/master/code/transcriptomically_similar_analysis) Mutational/Pathway analysis using transcriptomically similar samples (shared genes/recurrent alterations plots) |      Y      |    Y    |  Y  |     Y     |
|  [Immune analysis](https://github.com/d3b-center/OMPARE/tree/master/code/immune_analysis)                                                                                                                                    |      Y      |    N    |  N  |     N     |
|  [Survival analysis](https://github.com/d3b-center/OMPARE/tree/master/code/survival_analysis)                                                                                                                                  |      Y      |    N    |  N  |     N     |
|  [Genomic landscape plots](https://github.com/d3b-center/OMPARE/tree/master/code/genomic_landscape_plots) (circos, network)                                                                                                          |      Y      |    Y    |  Y  |     Y     |
|  [Oncogrid analysis](https://github.com/d3b-center/OMPARE/tree/master/code/oncogrid_analysis) (adapted from PNOC003)                                                                                                           |      Y      |    Y    |  Y  |     Y     |
|  [Drug recommendations](https://github.com/d3b-center/OMPARE/tree/master/code/drug_recommendations)                                                                                                                               |      Y      |    N    |  N  |     N     |
|  [Drug synergy](https://github.com/d3b-center/OMPARE/tree/master/code/drug_synergy)                                                                                                                                       |      Y      |    N    |  N  |     N     |

