# cohort3a + pnoc008 fpkm stranded
Rscript code/analyses/collapse-rnaseq/01-collapse-matrices.R \
--mat data/reference/hgg-dmg-integration/merged_ngs_files/rna-seq/cohort3a_wPNOC008-gene-expression-rsem-fpkm.stranded.rds \
--gene_sym FALSE \
--outfile cohort3a_wPNOC008-gene-expression-rsem-fpkm.stranded.collapsed.rds 

# cohort3a + pnoc008 tpm stranded
Rscript code/analyses/collapse-rnaseq/01-collapse-matrices.R \
--mat data/reference/hgg-dmg-integration/merged_ngs_files/rna-seq/cohort3a_wPNOC008-gene-expression-rsem-tpm.stranded.rds \
--gene_sym FALSE \
--outfile cohort3a_wPNOC008-gene-expression-rsem-tpm.stranded.collapsed.rds

# pbta + pnoc008 fpkm stranded
Rscript code/analyses/collapse-rnaseq/01-collapse-matrices.R \
--mat data/reference/hgg-dmg-integration/merged_ngs_files/rna-seq/pbta_wPNOC008-gene-expression-rsem-fpkm.stranded.rds \
--gene_sym FALSE \
--outfile pbta_wPNOC008-gene-expression-rsem-fpkm.stranded.collapsed.rds

# pbta + pnoc008 tpm stranded
Rscript code/analyses/collapse-rnaseq/01-collapse-matrices.R \
--mat data/reference/hgg-dmg-integration/merged_ngs_files/rna-seq/pbta_wPNOC008-gene-expression-rsem-tpm.stranded.rds \
--gene_sym FALSE \
--outfile pbta_wPNOC008-gene-expression-rsem-tpm.stranded.collapsed.rds
