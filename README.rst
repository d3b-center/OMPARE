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

	# NOTE: ggnetwork v0.5.1 is required

3. Download reference data:
   
.. code-block:: bash

	# get reference data from s3
	aws s3 --profile saml sync s3://d3b-bix-dev-data-bucket/PNOC008/reference /path-to/OMPARE/data/reference

4. Download patient-specific files from `data delivery project <https://cavatica.sbgenomics.com/u/cavatica/sd-8y99qzjj>`_:

* Copy Number: 

  * ``{uuid}.CNVs.p.value.txt``
  * ``{uuid}.controlfreec.info.txt``
  * ``{uuid}.controlfreec.ratio.txt``
  * ``{uuid}.gainloss.txt``

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

5. Download the following clinical and sample information files and add data manually to ``data/reference/manifest/manifest.xlsx`` 
   
* Files from `Kids First DRC <https://data-tracker.kidsfirstdrc.org/study/SD_8Y99QZJJ/documents>`_

  * PNOC008 Clinical Manifest
  * PNOC008 Sample Manifest

* Files from ADAPT (this part has been automated using ``code/update_pbta.R``): 
  
  * pbta-histologies-base-adapt.tsv

Detailed instructions are given in `d3b-analysis-toolkit <https://github.com/d3b-center/d3b-analysis-toolkit>`_

.. code-block:: bash

	cd /path-to/d3b-analysis-toolkit/scripts
	source .envrc
	python select-all-pbta-histologies.py -o /path-to/OMPARE/data/reference/pbta/pbta-histologies-base-adapt.tsv 

Scripts
=======

Master script
-------------

**run_OMPARE.R**: Master script that runs the following scripts:
   
1. **code/create_project_dir.R**: create project directory and organize files.
2. **code/create_clinfile.R**: create clinical file for patient of interest.
3. **code/update_pbta.R**: pull pbta histologies data from datawarehouse
4. **code/patient_level_analyses/pnoc_format.R**: update PNOC008 data matrices (cnv, mutations, fusions, expression) with each new patient.
5. **code/patient_level_analyses/gsea_enrichment.R**: update GSEA enrichment outputs with each new patient.
6. **code/patient_level_analyses/enrichment_output.R**: generate genes and pathway enrichment output for each new patient.
7. **OMPARE.Rmd**: run html reports
8. Using ``aws s3 sync``, sync back updated reference folder to ``s3://d3b-bix-dev-data-bucket/PNOC008/reference``
9. **upload_reports.R**: upload reports and output folders to PNOC008 data delivery project on cavatica.

.. code-block:: bash
	
	Options:
	--patient=PATIENT
		Patient identifier (PNOC008-22, C3342894...)

	--sourcedir=SOURCEDIR
		Source directory with all files (usually downloads folder on local machine)

	--clin_file=CLIN_FILE
		Manifest file (.xlsx)

	--update_pbta=UPDATE_PBTA
		Update PBTA adapt file (TRUE or FALSE)

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
	--clin_file data/reference/manifest/pnoc008_manifest.xlsx \
	--update_pbta FALSE \
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

**code/create_clinfile.R**: this script creates clinical file for patient of interest and stores under ``results/PNOC008-patient_num/clinical/``.

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
	--sheet /path-to/OMPARE/data/reference/manifest/pnoc008_manifest.xlsx \
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
	│   ├── {uuid}.controlfreec.CNVs.p.value.txt
	│   ├── {uuid}.controlfreec.info.txt
	│   ├── {uuid}.controlfreec.ratio.txt
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

**code/patient_level_analyses/pnoc_format.R**: this script updates the 008 patient matrices (cnv, mutations, fusions, expression) by adding current patient of interest
   
.. code-block:: bash

	Rscript code/patient_level_analyses/pnoc_format.R

	# Running the script will update the following files:
	data/reference/pnoc008
	├── pnoc008_clinical.rds
	├── pnoc008_cnv_filtered.rds
	├── pnoc008_consensus_mutation_filtered.rds
	├── pnoc008_counts_matrix.rds
	├── pnoc008_fpkm_matrix.rds
	├── pnoc008_fusions_filtered.rds
	├── pnoc008_tmb_scores.rds
	├── pnoc008_tpm_matrix.rds
	└── pnoc008_vs_gtex_brain_degs.rds


Update GSEA enrichment:
-----------------------

**code/patient_level_analyses/gsea_enrichment.R**: this script will update GSEA enrichment output with each new patient data.
   
.. code-block:: bash

	Rscript code/patient_level_analyses/gsea_enrichment.R --help

	Options:
	--patient=PATIENT
		Patient identifier for e.g. PNOC008-1, PNOC008-10 etc

	# Example for patient PNOC008-40
	Rscript code/patient_level_analyses/gsea_enrichment.R \
	--patient PNOC008-40

	# Running the script will update the following files:

	# reactome msigdb
	data/reference/gsea
	├── pbta_vs_gtex_brain.rds
	├── pbta_vs_pbta.rds
	├── pbta_vs_pbta_hgg.rds
	├── pnoc008_vs_gtex_brain.rds
	├── pnoc008_vs_pbta.rds
	├── pnoc008_vs_pbta_hgg.rds
	├── pnoc008_vs_tcga_gbm.rds
	├── tcga_gbm_vs_gtex_brain.rds
	└── tcga_gbm_vs_tcga_gbm.rds

	# dsigdb
	data/reference/dsigdb
	├── pnoc008_vs_gtex_brain.rds
	├── pnoc008_vs_pbta.rds
	└── pnoc008_vs_pbta_hgg.rds


Excel file with differential results:
-------------------------------------

**code/patient_level_analyses/enrichment_output.R**: this script will create an text file summaries containing up/down pathways and genes of patient of interest vs ``GTEx Brain``, ``PBTA HGG`` and ``PBTA all histologies``:

.. code-block:: bash

	Rscript code/patient_level_analyses/enrichment_output.R --help

	Options:
	--input=INPUT
		Directory e.g. data/PNOC008-04

	--output=OUTPUT
		output excel filename i.e. PNOC008-04_summary

	--type=TYPE
		text or excel

	# Example for patient PNOC008-40
	Rscript code/enrichment_output.R \
	--input /path-to/OMPARE/results/PNOC008-40 \
	--output PNOC008-40_summary \
	--type text

HTML reports:
-------------

Generate markdown report:

.. code-block:: bash

	# topDir is the project directory of current patient
	# fusion_method is the fusion method. Allowed values: star, arriba, both or not specified. (Optional) 
	# set_title is the title for the report. (Optional)
	# snv_pattern is one of the six values for simple variants: lancet, mutect2, strelka2, vardict, consensus, all (all four callers together)
	for(i in 1:length(callers)) {
    	output_dir <- file.path(topDir, 'Reports')
    	output_file <- paste0(patient, '_', callers[i], '.html')
    	input_file <- file.path(root_dir, 'OMPARE.Rmd')
    	rmarkdown::render(input = input_file,
    		params = list(topDir = topDir,
    			fusion_method = 'arriba',
        		set_title = set_title,
        		snv_caller = callers[i]), 
        	output_dir = output_dir, 
			intermediates_dir = output_dir,
			output_file = output_file)
	}


After running the reports, the project folder will have all output files with plots and tables under ``output`` and all html reports under ``reports``:

.. code-block:: bash

	results/PNOC008-40
	├── CEMITools
	│   ├── beta_r2.pdf
	│   ├── clustered_samples.rds
	│   ├── diagnostics.html
	│   ├── enrichment_es.tsv
	│   ├── enrichment_nes.tsv
	│   ├── enrichment_padj.tsv
	│   ├── expected_counts_corrected.rds
	│   ├── gsea.pdf
	│   ├── hist.pdf
	│   ├── hubs.rds
	│   ├── interaction.pdf
	│   ├── interactions.tsv
	│   ├── mean_k.pdf
	│   ├── mean_var.pdf
	│   ├── module.tsv
	│   ├── modules_genes.gmt
	│   ├── ora.pdf
	│   ├── ora.tsv
	│   ├── parameters.tsv
	│   ├── profile.pdf
	│   ├── qq.pdf
	│   ├── report.html
	│   ├── sample_tree.pdf
	│   ├── selected_genes.txt
	│   ├── summary.rds
	│   ├── summary_eigengene.tsv
	│   ├── summary_mean.tsv
	│   └── summary_median.tsv
	├── clinical
	│   └── patient_report.txt
	├── copy-number-variations
	│   ├── 106762e7-e100-405b-9ae9-bb80a186cdf9.controlfreec.CNVs.p.value.txt
	│   ├── 106762e7-e100-405b-9ae9-bb80a186cdf9.controlfreec.ratio.txt
	│   └── 106762e7-e100-405b-9ae9-bb80a186cdf9.gainloss.txt
	├── gene-expressions
	│   └── 806668be-e3a2-4ea3-90fb-f67eba78c7b3.rsem.genes.results.gz
	├── gene-fusions
	│   ├── 806668be-e3a2-4ea3-90fb-f67eba78c7b3.STAR.fusion_predictions.abridged.coding_effect.tsv
	│   └── 806668be-e3a2-4ea3-90fb-f67eba78c7b3.arriba.fusions.tsv
	├── output
	│   ├── PNOC008-40_summary_DE_Genes_Down.txt
	│   ├── PNOC008-40_summary_DE_Genes_Up.txt
	│   ├── PNOC008-40_summary_Pathways_Down.txt
	│   ├── PNOC008-40_summary_Pathways_Up.txt
	│   ├── circos_plot.png
	│   ├── cnv_plot.png
	│   ├── complexheatmap_cgs.png
	│   ├── complexheatmap_oncogrid.png
	│   ├── complexheatmap_phgg.png
	│   ├── consensus_mpf_output.txt
	│   ├── diffexpr_genes_barplot_output.rds
	│   ├── diffreg_pathways_barplot_output.rds
	│   ├── dim_reduction_plot_adult.rds
	│   ├── dim_reduction_plot_pediatric.rds
	│   ├── drug_dge_density_plots
	│   ├── drug_pathways_barplot.rds
	│   ├── dsigdb_de_genes_down.txt
	│   ├── dsigdb_de_genes_up.txt
	│   ├── dsigdb_pathways_down.txt
	│   ├── dsigdb_pathways_up.txt
	│   ├── filtered_germline_vars.rds
	│   ├── kaplan_meier_adult.rds
	│   ├── kaplan_meier_pediatric.rds
	│   ├── mutational_analysis_adult.rds
	│   ├── mutational_analysis_pediatric.rds
	│   ├── oncokb_cnv.txt
	│   ├── oncokb_cnv_annotated.txt
	│   ├── oncokb_consensus_annotated.txt
	│   ├── oncokb_fusion.txt
	│   ├── oncokb_fusion_annotated.txt
	│   ├── oncokb_lancet_annotated.txt
	│   ├── oncokb_merged_all_annotated.txt
	│   ├── oncokb_merged_all_annotated_actgenes.txt
	│   ├── oncokb_merged_consensus_annotated.txt
	│   ├── oncokb_merged_consensus_annotated_actgenes.txt
	│   ├── oncokb_merged_lancet_annotated.txt
	│   ├── oncokb_merged_lancet_annotated_actgenes.txt
	│   ├── oncokb_merged_mutect2_annotated.txt
	│   ├── oncokb_merged_mutect2_annotated_actgenes.txt
	│   ├── oncokb_merged_strelka2_annotated.txt
	│   ├── oncokb_merged_strelka2_annotated_actgenes.txt
	│   ├── oncokb_merged_vardict_annotated.txt
	│   ├── oncokb_merged_vardict_annotated_actgenes.txt
	│   ├── oncokb_mutect2_annotated.txt
	│   ├── oncokb_strelka2_annotated.txt
	│   ├── oncokb_vardict_annotated.txt
	│   ├── ora_plots.png
	│   ├── pathway_analysis_adult.rds
	│   ├── pathway_analysis_pediatric.rds
	│   ├── pbta_pnoc008_umap_output.rds
	│   ├── rnaseq_analysis_output.rds
	│   ├── ssgsea_scores_pediatric.rds
	│   ├── tcga_pnoc008_umap_output.rds
	│   ├── tis_scores.rds
	│   ├── tmb_profile_output.rds
	│   ├── transciptomically_similar_adult.rds
	│   ├── transciptomically_similar_pediatric.rds
	│   ├── transcriptome_drug_rec.rds
	│   └── tumor_signature_output.rds
	├── reports
	│   ├── PNOC008-40_all.html
	│   ├── PNOC008-40_consensus.html
	│   ├── PNOC008-40_lancet.html
	│   ├── PNOC008-40_mutect2.html
	│   ├── PNOC008-40_strelka2.html
	│   └── PNOC008-40_vardict.html
	└── simple-variants
	    ├── 106762e7-e100-405b-9ae9-bb80a186cdf9.lancet_somatic.vep.maf
	    ├── 106762e7-e100-405b-9ae9-bb80a186cdf9.mutect2_somatic.vep.maf
	    ├── 106762e7-e100-405b-9ae9-bb80a186cdf9.strelka2_somatic.vep.maf
	    ├── 106762e7-e100-405b-9ae9-bb80a186cdf9.vardict_somatic.vep.maf
	    ├── c185fc36-97d9-433d-9ea9-25a608b2f660.gatk.PASS.vcf.gz.hg38_multianno.txt.gz
	    └── e9248ac8-79e5-41e7-a97d-3ccd9c406074.consensus_somatic.vep.maf

	9 directories, 105 files


Upload to data-delivery project
-------------------------------

**upload_reports.R**: this script uploads the files under ``reports``, ``output`` and ``CEMITools`` folder to the data delivery project folder on cavatica. 

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

