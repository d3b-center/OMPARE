## Script Author

Run Jin

## Calculate synergy scores for combination of drug recommendations
This module is written to implement the synergy from gene expression and network mining (SynGeNet) method described in this paper(https://www.nature.com/articles/s41540-019-0085-4) with custom modifications.

Briefly, we leveraged the outputs from CEMiTools and signatureSearch (which are already incorporated in our OMPARE workflow) and calculate synergy score as followed:
1. From `clustered_samples.rds` and `pnoc008_clinical.rds`, extract the cluster number that is assigned to our sample of interest 
2. From `enrichment_nes.tsv`, find **all** the module numbers that are positively correlated with the cluster of our sample of interest.
3. Generate subnetwork associated with the modules found above by subsetting `interactions.tsv` to contain only interactions in that particular module - result saved as `reults/subnetwork_genes.tsv`
4. Use bioconductor package `drugTargetInteractions` and `biomaRt` to find all FDA approved drugs that target genes in our subnetwork, annotate with drug actions, MOA etc. - results saved as  `results/drug_gene_map/subnetwork_gene_drug_map.tsv`.
5. Additionally, drugs listed in *.qSig.txt from signatureSearch (which are drugs whose gene expression signatures are _reversely_ correlated with gene expression signature of our patients of interest) were intersected with each module listed the `results/drug_gene_map/subnetwork_gene_drug_map.tsv`
Results save as `results/drug_gene_map/*_qSig_subnetwork_drug_gene_map.tsv`.
**NOTE**: for each patient of interest, 3 gene expression signature are available: DGE compared to GTEx brain normal, DGE compared to all PBTA samples and DGE compared to PBTA HGG samples. 

After drug gene mapped results were saved, these results were used to calculate drug synergy score as followed:
1. Firstly, the subnetwork and mapped qSig files were filtered to one particular module 
2. Module-subsetted subnetwork was graphed using `graph.edgelist` function from `igraph` to generate a graph object, which will be used as input for calculating synergy score.
3. All the drugs listed in module-subsetted `*_qSig_subnetwork_drug_gene_map.tsv` file will be given a weighted score between 1-2 based on their WTCS score - the lower the score (more negative), the higher the rank.
4. For each drug, we then subset their targets to only targets that are in the module-subsetted subnetwork, and then we calculate the centrality score of that drug based on the `closeness` and `betweenness` of its targets in the module-subsetted subnetwork.
5. Synergy score of a drug combination drugA and drugB is calculated as: weighted score of drugA (from step 2) * weighted score of drugB (from step 2) * centrality score of all targets from both drugs in the subnetwork. 
6. Same calculation (step 1-5) was done for all other modules in the subnetwork. 
7. Drug combinations along with calculated synergy score are output and saved as: `results/synergy_score/*_qSig_synergy_score.tsv`.

### Usage

All data files required are:
  - For the patient of interest (download from Cavatica delivery project):
    - `interactions.tsv`
    - `enrichment_nes.tsv`
    - `GTExBrain_qSig_output.txt`
    - `PBTA_ALL_qSig_output.txt`
    - `PBTA_HGG_qSig_output.txt`
    - `clustered_samples.rds`
  
All reference files required are:
  - References CHEMBL database (download from https://chembl.gitbook.io/chembl-interface-documentation/downloads):
    - `chembl_29.db`
  - Clinical files for all PNOC008 samples (download from S3 bucket s3://d3b-bix-dev-data-bucket/PNOC008/reference/pnoc008/)
    - `pnoc008_clinical.rds`
    
After all the relevant files are downloaded and stored in either `data` or `references` folder of the main repository, the module can be run with:
  
  ```
bash run_synergy.sh
```

###### Contents

`01-subnetwork_qSig_gene_drug_map.R` generates drug gene map for subnetwork of interest (subnetwork that is **most** positively correalted with the cluster assigned to our sample). 
Additionally, drugs identified by the *qSig output of signatureSearch were subset to drugs that have targets in the subnetwork of interest and annotated with MOA, drug actions etc by `drugTargetInteractions` package.

Input:
  - `../../data/interactions.tsv`
  - `../../data/enrichment_nes.tsv`
  - `../../data/clustered_samples.rds`
  - `../../data/GTExBrain_qSig_output.txtv`
  - `../../data/PBTA_ALL_qSig_output.txt`
  - `../../data/PBTA_HGG_qSig_output.txt`
  - `../../references/chembl_29_sqlite/chembl_29.db`
  - `../../references/pnoc008_clinical.rds`

Output:
  - `results/subnetwork_genes.tsv`
  - `results/drug_gene_map/subnetwork_gene_drug_map.tsv`
  - `results/drug_gene_map/gtex_qSig_subnetwork_drug_gene_map.tsv`
  - `results/drug_gene_map/pbta_qSig_subnetwork_drug_gene_map.tsv`
  - `results/drug_gene_map/pbta_hgg_qSig_subnetwork_drug_gene_map.tsv`
  
`02-drug_synergy_score_calc.R` calculates the synergy score using the method described above.

Input:
  - `results/subnetwork_genes.tsv`
  - `results/drug_gene_map/gtex_qSig_subnetwork_drug_gene_map.tsv`
  - `results/drug_gene_map/pbta_qSig_subnetwork_drug_gene_map.tsv`
  - `results/drug_gene_map/pbta_hgg_qSig_subnetwork_drug_gene_map.tsv`

Output:
  - `results/synergy_score/gtex_qSig_synergy_score.tsv`
  - `results/synergy_score/pbta_qSig_synergy_score.tsv`
  - `results/synergy_score/pbta_hgg_qSig_synergy_score.tsv`
  - `results/synergy_score/combined_qSig_synergy_score.tsv`


