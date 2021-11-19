.. |date| date::

********************
Omics Patient Report
********************

:authors: Komal S. Rathi, Adam Kraya, Run Jin
:contact: RATHIK@chop.edu, KRAYAA@chop.edu, JINR@chop.edu
:organization: D3B, CHOP
:status: This is "work in progress"
:date: |date|

.. meta::
   :keywords: omics, report, flexboard, 2019
   :description: Omics Patient Report

Prerequisites
=============

1. Clone the OMPARE repository.

2. Install R packages:

.. code-block:: bash

	# install packages
	cd /path-to/OMPARE
	Rscript code/utils/install_pkgs.R

	# NOTE: ggnetwork v0.5.1 and DT.0.18 is required

3. Download reference data:
   
.. code-block:: bash

	# get reference data from s3
	aws s3 --profile saml sync s3://d3b-bix-dev-data-bucket/PNOC008/reference /path-to/OMPARE/data/

4. Download patient-specific files from `data delivery project <https://cavatica.sbgenomics.com/u/cavatica/sd-8y99qzjj>`_:

* Copy Number: 

  * ``{uuid}.controlfreec.info.txt`` (for purity and ploidy)
  * ``{uuid}.gainloss.txt``
  * ``{uuid}.diagram.pdf``

* Expression:

  * ``{uuid}.genes.results``

* Fusions: 

  * ``{uuid}.arriba.fusions.tsv``
  * ``{uuid}.STAR.fusion_predictions.abridged.coding_effect.tsv``

* Somatic Variants: 
 
  * ``{uuid}.{lancet, mutect2, strelka2, vardict}_somatic.norm.annot.protected.maf``
  * ``{uuid}.consensus_somatic.protected.maf``

* Germline Variants: 

  * ``{uuid}.gatk.PASS.vcf.gz.hg38_multianno.txt.gz``

5. Download the following clinical information files and add data manually to ``data/manifest/manifest.xlsx`` 
   
* Files from `Kids First DRC <https://data-tracker.kidsfirstdrc.org/study/SD_8Y99QZJJ/documents>`_

  * PNOC008 Clinical Manifest (needed to map ``Research ID`` to ADAPT ``cohort_participant_id``)

* Files from ADAPT (updated each morning and needed to get BS identifier and other information)

.. code-block:: bash

	aws s3 --profile saml cp s3://d3b-bix-dev-data-bucket/pbta-histologies-base-adapt.tsv data/pbta/

Note: None of these files have information on short_histology or broad_histology so currently it is being hard-coded ``HGAT`` and ``Diffuse astrocytic and oligodendroglial tumor``, respectively.


Scripts
=======

Master script
-------------

**run_OMPARE.R**: Master script that runs the following scripts:
   
1. **code/create_project_dir.R**: create project directory and organize files.
2. **code/create_clinfile.R**: create clinical file for patient of interest.
3. **code/update_pnoc008_matrices.R**: update PNOC008 data matrices (cnv, mutations, fusions, expression) with each new patient.
4. **OMPARE.Rmd**: run html reports
5. Using ``aws s3 sync``, sync back updated data folder to ``s3://d3b-bix-dev-data-bucket/PNOC008/reference``
6. **upload_reports.R**: upload reports and output folders to PNOC008 data delivery project on cavatica.

.. code-block:: bash
	
	Options:
	--patient=PATIENT
		Patient identifier (PNOC008-22, C3342894...)

	--source_dir=SOURCE_DIR
		Source directory with all files

	--clin_file=CLIN_FILE
		Manifest file (.xlsx)

	--sync_data=SYNC_DATA
		Sync reference data to s3 (TRUE or FALSE)

	--upload_reports=UPLOAD_REPORTS
		Upload reports to cavatica (TRUE or FALSE)

	--study=STUDY
		Study ID (PNOC008 or CBTN)

	# Example for patient PNOC008-40
	Rscript run_OMPARE.R \
	--patient PNOC008-40 \
	--sourcedir ~/Downloads/p40 \
	--clin_file data/manifest/pnoc008_manifest.xlsx \
	--sync_data TRUE \
	--upload_reports FALSE \
	--study PNOC008

Create project directory
------------------------

**code/create_project_dir.R**: this script creates and organizes input files under ``results``. Creates ``output`` folder to store all output for plots and tables reported and ``reports`` folder to store all html output.
   
.. code-block:: bash

	Rscript code/create_project_dir.R --help

	Options:
	--sourcedir=SOURCEDIR
		Source directory with all files

	--destdir=DESTDIR
		Destination directory.

	# Example for patient PNOC008-40
	Rscript code/create_project.R \
	--sourcedir ~/Downloads/p40 \
	--destdir /path-to/OMPARE/results/PNOC008-40

Create clinical file
--------------------

**code/create_clinfile.R**: this script creates clinical file for patient of interest and stores under ``results/PNOC008-XX/clinical/``.

.. code-block:: bash

	Rscript code/create_clinfile.R --help

	Options:
	--sheet=SHEET
		PNOC008 Manifest file (.xlsx)

	--dir=DIR
		Path to PNOC008 patient folder.

	--patient=PATIENT
		Patient identifier for PNOC008. e.g. PNOC008-1, PNOC008-10 etc

	# Example for patient PNOC008-40
	Rscript code/create_clinfile.R \
	--sheet /path-to/OMPARE/data/manifest/pnoc008_manifest.xlsx \
	--patient PNOC008-40 \
	--dir /path-to/OMPARE/results/PNOC008-40

NOTE: The above steps will create a directory structure for the patient of interest: 

.. code-block:: bash

	# Example for PNOC008-40
	.
	results/PNOC008-40
	├── clinical
	│   └── patient_report.txt
	├── copy-number-variations
	│   ├── {uuid}.controlfreec.info.txt
	│   ├── {uuid}.diagram.pdf	
	│   └── {uuid}.gainloss.txt
	├── gene-expressions
	│   └── {uuid}.rsem.genes.results.gz
	├── gene-fusions
	│   ├── {uuid}.STAR.fusion_predictions.abridged.coding_effect.tsv
	│   └── {uuid}.arriba.fusions.tsv
	├── output
	├── reports
	└── simple-variants
	    ├── {uuid}.lancet_somatic.norm.annot.protected.maf
	    ├── {uuid}.mutect2_somatic.norm.annot.protected.maf
	    ├── {uuid}.strelka2_somatic.norm.annot.protected.maf
	    ├── {uuid}.vardict_somatic.norm.annot.protected.maf
	    ├── {uuid}.consensus_somatic.protected.maf
	    └── {uuid}.gatk.PASS.vcf.gz.hg38_multianno.txt.gz


Update PNOC008 data matrices:
-----------------------------

**code/update_pnoc008_matrices.R**: this script updates the 008 patient matrices (cnv, mutations, fusions, expression) by adding current patient of interest
   
.. code-block:: bash

	Rscript code/update_pnoc008_matrices.R

	# Running the script will update the following files:
	data/pnoc008
	├── pnoc008_clinical.rds
	├── pnoc008_cnv_filtered.rds
	├── pnoc008_consensus_mutation_filtered.rds
	├── pnoc008_counts_matrix.rds
	├── pnoc008_fpkm_matrix.rds
	├── pnoc008_fusions_filtered.rds
	├── pnoc008_tmb_scores.rds
	├── pnoc008_tpm_matrix.rds
	└── pnoc008_vs_gtex_brain_degs.rds

HTML reports:
-------------

Generate markdown report:

.. code-block:: bash

	# patient_dir is the project directory of current patient
	# set_title is the title for the report. (Optional)
	# snv_pattern is one of the six values for simple variants: lancet, mutect2, strelka2, vardict, consensus, all (all four callers together)
	Rscript -e "rmarkdown::render(input = 'OMPARE.Rmd', 
	params = list(patient_dir = patient_dir,
			set_title = set_title,
			snv_caller = snv_caller), 
			output_dir = output_dir, 
			intermediates_dir = output_dir,
			output_file = output_file, clean = TRUE)"

After running the reports, the project folder will have all output files with plots and tables under ``output`` and all html reports under ``reports``:

.. code-block:: bash

	.
	├── drug_recommendations
	│   ├── CEMITools
	│   │   ├── beta_r2.pdf
	│   │   ├── clustered_samples.rds
	│   │   ├── diagnostics.html
	│   │   ├── enrichment_es.tsv
	│   │   ├── enrichment_nes.tsv
	│   │   ├── enrichment_padj.tsv
	│   │   ├── expected_counts_corrected.rds
	│   │   ├── gsea.pdf
	│   │   ├── hist.pdf
	│   │   ├── hubs.rds
	│   │   ├── interaction.pdf
	│   │   ├── interactions.tsv
	│   │   ├── mean_k.pdf
	│   │   ├── mean_var.pdf
	│   │   ├── module.tsv
	│   │   ├── modules_genes.gmt
	│   │   ├── ora.pdf
	│   │   ├── ora.tsv
	│   │   ├── parameters.tsv
	│   │   ├── profile.pdf
	│   │   ├── qq.pdf
	│   │   ├── report.html
	│   │   ├── sample_tree.pdf
	│   │   ├── selected_genes.txt
	│   │   ├── summary.rds
	│   │   ├── summary_eigengene.tsv
	│   │   ├── summary_mean.tsv
	│   │   ├── summary_median.tsv
	│   │   ├── umap_output.rds
	│   │   └── umap_top_20_neighbors_output.rds
	│   ├── GTExBrain_dsea_go_mf_output.html
	│   ├── GTExBrain_dsea_go_mf_output.pdf
	│   ├── GTExBrain_dsea_go_mf_output.txt
	│   ├── GTExBrain_dsea_go_mf_output_files
	│   ├── GTExBrain_qSig_output.txt
	│   ├── GTExBrain_tsea_reactome_output.txt
	│   ├── PBTA_ALL_dsea_go_mf_output.html
	│   ├── PBTA_ALL_dsea_go_mf_output.pdf
	│   ├── PBTA_ALL_dsea_go_mf_output.txt
	│   ├── PBTA_ALL_dsea_go_mf_output_files
	│   ├── PBTA_ALL_qSig_output.txt
	│   ├── PBTA_ALL_tsea_reactome_output.txt
	│   ├── PBTA_HGG_dsea_go_mf_output.html
	│   ├── PBTA_HGG_dsea_go_mf_output.pdf
	│   ├── PBTA_HGG_dsea_go_mf_output.txt
	│   ├── PBTA_HGG_dsea_go_mf_output_files
	│   ├── PBTA_HGG_qSig_output.txt
	│   ├── PBTA_HGG_tsea_reactome_output.txt
	│   ├── {patient_id}_CHEMBL_drug-gene.tsv
	│   ├── drug_dge_density_plots
	│   │   ├── {gene}_drug_dge_density_plots.png
	│   │   └── top_drug_dge_density_plots.pdf
	│   ├── drug_pathways_barplot.pdf
	│   ├── ora_plots.pdf
	│   └── transcriptome_drug_rec.rds
	├── drug_synergy
	│   ├── combined_qSig_synergy_score.tsv
	│   ├── combined_qSig_synergy_score_top10.pdf
	│   ├── gtex_qSig_subnetwork_drug_gene_map.tsv
	│   ├── gtex_qSig_synergy_score.tsv
	│   ├── pbta_hgg_qSig_subnetwork_drug_gene_map.tsv
	│   ├── pbta_hgg_qSig_synergy_score.tsv
	│   ├── pbta_qSig_subnetwork_drug_gene_map.tsv
	│   ├── pbta_qSig_synergy_score.tsv
	│   ├── subnetwork_gene_drug_map.tsv
	│   └── subnetwork_genes.tsv
	├── filtered_germline_vars.rds
	├── genomic_landscape_plots
	│   └── circos_plot.png
	├── immune_analysis
	│   ├── immune_scores_adult.pdf
	│   ├── immune_scores_adult.rds
	│   ├── immune_scores_pediatric.pdf
	│   ├── immune_scores_pediatric.rds
	│   ├── immune_scores_topcor_pediatric.pdf
	│   ├── immune_scores_topcor_pediatric.rds
	│   ├── tis_scores.pdf
	│   └── tis_scores.rds
	├── oncogrid_analysis
	│   └── complexheatmap_oncogrid.pdf
	├── oncokb_analysis
	│   ├── oncokb_cnv.txt
	│   ├── oncokb_cnv_annotated.txt
	│   ├── oncokb_fusion.txt
	│   ├── oncokb_fusion_annotated.txt
	│   ├── oncokb_{snv_caller}_annotated.txt
	│   ├── oncokb_merged_{snv_caller}_annotated.txt
	│   └── oncokb_merged_{snv_caller}_annotated_actgenes.txt
	├── rnaseq_analysis
	│   ├── {patient_id}_summary_DE_Genes_Down.txt
	│   ├── {patient_id}_summary_DE_Genes_Up.txt
	│   ├── {patient_id}_summary_Pathways_Down.txt
	│   ├── {patient_id}_summary_Pathways_Up.txt
	│   ├── diffexpr_genes_barplot_output.rds
	│   ├── diffreg_pathways_barplot_output.rds
	│   └── rnaseq_analysis_output.rds
	├── survival_analysis
	│   ├── kaplan_meier_adult.pdf
	│   └── kaplan_meier_pediatric.pdf
	├── tmb_analysis
	│   ├── consensus_mpf_output.txt
	│   ├── tmb_profile_output.rds
	│   └── tumor_signature_output.rds
	└── transcriptomically_similar_analysis
	    ├── dim_reduction_plot_adult.rds
	    ├── dim_reduction_plot_pediatric.rds
	    ├── lollipop_recurrent_adult.pdf
	    ├── lollipop_recurrent_pediatric.pdf
	    ├── lollipop_shared_adult.pdf
	    ├── lollipop_shared_pediatric.pdf
	    ├── mutational_analysis_adult.rds
	    ├── mutational_analysis_pediatric.rds
	    ├── mutational_cnv_recurrent_adult.pdf
	    ├── mutational_cnv_recurrent_pediatric.pdf
	    ├── mutational_cnv_shared_adult.pdf
	    ├── mutational_cnv_shared_pediatric.pdf
	    ├── mutational_recurrent_adult.pdf
	    ├── mutational_recurrent_pediatric.pdf
	    ├── mutational_shared_adult.pdf
	    ├── mutational_shared_pediatric.pdf
	    ├── pathway_analysis_adult.pdf
	    ├── pathway_analysis_adult.rds
	    ├── pathway_analysis_pediatric.pdf
	    ├── pathway_analysis_pediatric.rds
	    ├── pbta_hgat_pnoc008_nn_table.rds
	    ├── pbta_hgat_pnoc008_umap_output.rds
	    ├── pbta_pnoc008_nn_table.rds
	    ├── pbta_pnoc008_umap_output.rds
	    ├── ssgsea_scores_pediatric.pdf
	    ├── ssgsea_scores_pediatric.rds
	    ├── tcga_gbm_pnoc008_nn_table.rds
	    ├── tcga_pnoc008_umap_output.rds
	    ├── transciptomically_similar_adult.rds
	    └── transciptomically_similar_pediatric.rds



Upload to data-delivery project
-------------------------------

**upload_reports.R**: this script uploads the files under ``reports`` and ``output`` folders to the data delivery project folder on cavatica. 

.. code-block:: bash

	Rscript upload_reports.R --help

    Options:
	--patient=PATIENT
		Patient Identifier (PNOC008-22, etc...)

	--study=STUDY
		PNOC008 or CBTN

	# Example run for PNOC008-40
	Rscript upload_reports.R \
	--patient PNOC008-40 \
	--study 'PNOC008'

Dependencies on specific hgg-dmg versions
=========================================

These hgg-dmg files are ``20201202-data`` version dependent:

.. code-block:: bash

	hgg-dmg-integration
	└── 20201202-data
	    ├── CC_based_heatmap_Distance_euclidean_finalLinkage_average_clusterAlg_KM_expct_counts_VST_cluster_and_annotation.tsv
	    ├── pbta-hgat-dx-prog-pm-gene-counts-rsem-expected_count-uncorrected.rds
	    └── pbta-histologies.tsv

