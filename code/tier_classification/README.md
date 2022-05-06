## Tier and Hotspot Annotation for Samples Used in OMPARE

This module contains R script that added `Tier: xx` and `Cancer Hotspot` to the tables in P1:
`key_clinical_findings_output` and `all_findings_output`.

### Usage
The module can be run with the R script in the module:
  
```
Rscript run_tier_classification.R
```

The module annotate tier based on Oncokb annotations using tier_classification based on the following rules: 

Input file `oncokb_consensus_annotated.txt`

1. Tier 1 variants: 
  a. Levels of evidence to include among Tier 1 would be: in the `HIGHEST_LEVEL` field, those with levels 1-3 or R1; OR in the `Highest DX Level`, levels 1-2; OR in the `Highest PX Level`, levels 1-2; OR annotated as `Tier 1` in COSMIC resistance variant database; OR designated as `pathogenic` in `CLIN_SIG`.
  b. Mutation type: from MAF `Variant_Classification` field - missense, nonsense, indel (frameshift/non-frameshift), splice site, splice region.
  c. Variant allele frequency > 0.05: Calculated as `VAF = t_alt_count/(t_alt_count+t_ref_count)`
  d. Population database: from MAF - `gnomad_AF` < 0.01
  e. From MAF: `Existing_variation`: Must contain `COSM` identifier (found in COSMIC) 
  g. Pathways: Known TSG or Oncogene given OMPARE knowledgebase (`cancer_gene_list.rds`). 
  h. Publications: from `oncokb_consensus_annotated.txt`, `CITATIONS` is non-empty

2. Tier 2 variants: 
  a. Levels of evidence to include among Tier 2 would be: - in the `HIGHEST_LEVEL` field, those variants only with level 4 evidence OR R2; OR in the `Highest DX Level`, level 3; OR in the `Highest PX Level`, level 3; OR present in COSMIC resistance variant database without `Tier 1` annotated; OR designated as `likely_pathogenic` in `CLIN_SIG`.
(Note: a given variant may have multiple Levels of evidence in different indications, so long as the levels don't exceed 3 for DX or PX, nor 4 for `Level`, it would remain a Tier2 variant. If a variant has a level higher in any single one of these categories, it is a Tier 1 variant). 
  b. Mutation type: from MAF `Variant_Classification` field - missense, nonsense, indel (frameshift/non-frameshift), splice site, splice region.
  c. Variant allele frequency > 0.05: Calculated as `VAF = t_alt_count/(t_alt_count+t_ref_count)`
  d. Population database: from MAF - `gnomad_AF` < 0.01
  e. From MAF: `Existing_variation`: Must contain `COSM` identifier (found in COSMIC) 
  g. Pathways: Known TSG or Oncogene 
  h. Publications: from `oncokb_consensus_annotated.txt`, `CITATIONS` is non-empty

3. Tier 3 variants (VUS):
  a. No levels annotation in `oncokb_consensus_annotated.txt`
  b. Mutation type: from MAF `Variant_Classification` field - missense, in-frame insertions and deletions
  c. Variant allele frequency > 0.05: Calculated as `VAF = t_alt_count/(t_alt_count+t_ref_count)`
  d. Population database: from MAF - `gnomad_AF` < 0.01
  e. From MAF: `Existing_variation`: May or may not contain `COSM` identifier (found in COSMIC). 
  g. Pathways: Gene may or may not map be included in TSG or oncogene list
  h. Publications: from `oncokb_consensus_annotated.txt`, `CITATIONS` may or may not be empty

4. Tier 4 variants (Benign): 
  a. No levels annotation in `oncokb_consensus_annotated.txt`
  b. Mutation type: from MAF `Variant_Classification` field - missense
  c. Variant allele frequency either 40%-60% or >90% (non-mosaic): Calculated as `VAF = t_alt_count/(t_alt_count+t_ref_count)`
  d. Population database: from MAF - `gnomad_AF` >= 0.01
  e. From MAF:  `SIFT` and `Polyphen` do not contain `deleterious` or `damaging`, respectively
  f. Pathways: Gene is not a known TSG or oncogene
  g. Publications: from `oncokb_consensus_annotated.txt`, `CITATIONS` is empty
  
`hotspot_database_2017_indel.tsv` and `hotspot_database_2017_snv.tsv` derived from curated MSKCC data (downloaded from https://www.cancerhotspots.org/#/download), v2 version, using [this script](https://github.com/runjin326/CHOP_miscellaneous/blob/main/hotspot_prepare/hotspot_prep.R)

To see whether SNV hotspot were present, the `key_clinical_findings_output` and `all_findings_output` were first filtered to contain only `Missense_Mutation`, `Nonsense_Mutation`, `Splic_Variant` and `Splice_Region`. 
For every entry in the hotspot file, the gene symbol and AA position were used to query the `key_clinical_findings_output` and `all_findings_output` table to see whether they are present.
If the exact variant is present, then `Cancer Hotspot` annotation will be added to the `Variant_Properties` column;
If the exact variant is not present but the AA position matches with hotspot AA position, then then `Cancer Hotspot Location` annotation will be added to the `Variant_Properties` column.

To see whether indel hotspot were present, the `key_clinical_findings_output` and `all_findings_output` were first filtered to contain only `Frame_Shift_Del`, `Frame_Shift_Ins`m `In_Frame_Del` and `In_Frame_Ins`. 
For every entry in the indel hotspot file, the gene symbol and `HGVSp_Short` were used to query the `key_clinical_findings_output` and `all_findings_output` table to see whether they are present - if yes, then `Cancer Hotspot` annotation will be added to the `Variant_Properties` column.

This module then add CIVIC curated annotation to the `all_findings_output` table using the script `civic_tier_anno.R` with the following rules:

1. Annotate mutations
1) First annotate specific SNV using `c_dot_civic.tsv` and `snv_remaining_civic.tsv`
2) And then match by domain or exon and annotate any mutation within the domain/exon if they are in the CIVIC domain or exon files `mutation_exon_civic_annotated.tsv` and `mutation_domain_civic_annotated.tsv`
3) insertions or deletions that are within exons described by CIVIC will be annotated as well (files used are `del_exon_civic_annotated.tsv` and `ins_exon_civic_annotated.tsv`)
4) lastly, any remaining entries with mutations in genes listed in `mutation_general_civic.tsv` will be annotated 

2. Annotate fusions using `fusion_column_only.tsv` and `fusion_word_only.tsv`
1) `fusion_column_only.tsv` is used to annotate specific fusions with both gene partners
2) then if any fusion has a fusion partner in the file `fusion_word_only.tsv`, the entry will be annotated

3. Annotate expression results using `expression_all_civic.tsv`
4. Annotate gain and loss from CNV results using `amplification_only_civic.tsv` and `deletion_general_civic.tsv`

5. Annotate all the rest using `alteration_civic.tsv`
For any of the entries that were not annotated by the previous 4 steps, they would be annotated if the gene is annotated in `alteration_civic.tsv`.

Based on the civic evidence level, add additional tier or upgrade tier based on the following rules 
1. Tier I: designated as A or B in `evidence_level_civic`
2. Tier II: designated as C or D in `evidence_level_civic`
3. Tier III: designated as C or D in `evidence_level_civic`
