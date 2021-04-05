.. |date| date::

********************
Omics Patient Report
********************

:authors: Komal S Rathi, Pichai Raman
:contact: rathik@email.chop.edu
:organization: D3B, CHOP
:status: This is "work in progress"
:date: |date|

.. meta::
   :keywords: omics, report, flexboard, 2019
   :description: Omics Patient Report

Prerequisites
=============

1. Clone the repository.

2. Install R Packages:

.. code-block:: bash

	# install packages
	Rscript code/utils/install_pkgs.R

	# NOTE: ggnetwork v0.5.1 is required

3. Download reference Data:
   
.. code-block:: bash

	# get reference data from s3
	aws s3 --profile saml s3://d3b-bix-dev-data-bucket/PNOC008/Reference /path/to/OMPARE/data/reference

4. Download input files from data delivery project:

* Copy Number: ``*.CNVs.p.value.txt``
* Copy Number: ``*.controlfreec.ratio.txt``
* Copy Numberr: ``*.gainloss.txt``
* Expression: ``*.genes.results``
* Fusions: ``*.arriba.fusions.tsv``
* Fusions: ``*.STAR.fusion_predictions.abridged.coding_effect.tsv``
* Somatic Variants: ``*_somatic.vep.maf``
* Germline Variants: ``*.gatk.PASS.vcf.gz.hg38_multianno.txt.gz``

Scripts
=======

Master script
-------------

**run_OMPARE.R**: Master script that runs the following scripts:
   
1. **code/create_project_dir.R**: creating project directory and organize files.
2. **code/create_clinfile.R**: creating clinical file for patint of interest.
3. **code/patient_level_analyses/pnoc_format.R**: updating PNOC008 data matrices (cnv, mutations, fusions, expression) with each new patient.
4. **code/patient_level_analyses/gsea_enrichment.R**: updating GSEA enrichment outputs with each new patient.
5. **code/patient_level_analyses/enrichment_output.R**: generating genes and pathway enrichment output for each new patient.
6. **OMPARE.Rmd**: running html reports

.. code-block:: bash
	
	Rscript run_OMPARE.R --help

	Options:
	-p PATIENT, --patient=PATIENT
		Patient Number (1, 2...)

	-s SOURCEDIR, --sourcedir=SOURCEDIR
		Source directory with all files

	-c CLIN_FILE, --clin_file=CLIN_FILE
		PNOC008 Manifest file (.xlsx)

	# Example run for PNOC008-21
	Rscript run_OMPARE.R \
	--patient 21 \
	--clin_file data/reference/Manifest/PNOC008_Manifest.xlsx \
	--sourcedir /path/to/downloaded_files_from_cavatica


Create project directory
------------------------

**code/create_project_dir.R**: this script creates and organizes input files under ``results``. Creates ``output`` folder to store all output for plots and tables reported and ``reports`` folder to store all html output.
   
.. code-block:: bash

	Rscript code/create_project_dir.R --help

	Options:
		-s SOURCEDIR, --sourcedir=SOURCEDIR
			Source directory containing all files from data delivery project

		-d DESTDIR, --destdir=DESTDIR
			Destination directory. Should be /path/to/OMPARE/results/PNOC008-21/ for Patient 21

		-h, --help
			Show this help message and exit

	# Example for Patient PNOC008-21
	Rscript code/create_project.R \
	--sourcedir /path/to/source/PNOC008-21-cavatica-files \
	--destdir /path/to/OMPARE/results/PNOC008-21/

Create clinical file
--------------------

**code/create_clinfile.R**: this script creates clinical file for patient of interest and stores under ``results/PNOC008-patient_num/clinical/``.

.. code-block:: bash

	Rscript code/create_clinfile.R --help

	Options:
		-s SHEET, --sheet=SHEET
			PNOC008 Manifest file (.xlsx)

		-d DIR, --dir=DIR
			Path to PNOC008 patient folder.

		-p PATIENT, --patient=PATIENT
			Patient identifier for PNOC008. e.g. PNOC008-1, PNOC008-10 etc

	# Example for Patient PNOC008-21
	Rscript code/create_clinfile.R \
	--sheet data/reference/Manifest/PNOC008_Manifest.xlsx \
	--patient PNOC008-21 \
	--dir /path/to/OMPARE/results/PNOC008-21

NOTE: The above steps will create a directory structure for the patient of interest: 

.. code-block:: bash

	# Example for PNOC008-21
	.
	results/PNOC008-21
	├── clinical
	│   └── patient_report.txt
	├── copy-number-variations
	│   ├── uuid.controlfreec.CNVs.p.value.txt
	│   └── uuid.controlfreec.ratio.txt
	├── gene-expressions
	│   └── uuid.rsem.genes.results.gz
	├── gene-fusions
	│   ├── uuid.STAR.fusion_predictions.abridged.coding_effect.tsv
	│   └── uuid.arriba.fusions.tsv
	├── output
	├── reports
	└── simple-variants
	    ├── uuid.lancet_somatic.vep.maf
	    ├── uuid.mutect2_somatic.vep.maf
	    ├── uuid.strelka2_somatic.vep.maf
	    ├── uuid.vardict_somatic.vep.maf
	    ├── uuid.consensus_somatic.vep.maf
	    └── uuid.gatk.PASS.vcf.gz.hg38_multianno.txt.gz

Update GSEA enrichment:
-----------------------

**code/patient_level_analyses/gsea_enrichment.R**: this script will update GSEA enrichment output with each new patient data.
   
.. code-block:: bash

	Rscript code/patient_level_analyses/gsea_enrichment.R

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

Excel file with differential results:
-------------------------------------

**code/patient_level_analyses/enrichment_output.R**: this script will create an text file summaries containing up/down pathways and genes of patient of interest vs ``GTEx Brain``, ``PBTA HGG`` and ``PBTA all histologies``:

.. code-block:: bash

	Rscript code/patient_level_analyses/enrichment_output.R --help

	Options:
		-i INPUT, --input=INPUT
			Directory e.g. data/PNOC008-04

		-o OUTPUT, --output=OUTPUT
			output excel filename i.e. PNOC008-04_summary

		-t TYPE, --type=TYPE
			text or excel

	# Example for Patient PNOC008-21
	Rscript code/enrichment_output.R \
	--input /path/to/OMPARE/results/PNOC008-21 \
	--output PNOC008-21_summary \
	--type text

HTML reports:
-------------

8. Generate markdown report:

.. code-block:: bash

	# topDir is the project directory of current patient
	# fusion_method is the fusion method. Allowed values: star, arriba, both or not specified. (Optional) 
	# set_title is the title for the report. (Optional)
	# snv_pattern is one of the six values for simple variants: lancet, mutect2, strelka2, vardict, consensus, all (all four callers together)
	# tmb (Tumor mutational burden) is set to 77.46.
	for(i in 1:length(callers)) {
    	output_dir <- file.path(topDir, 'Reports')
    	output_file <- paste0(patient, '_', callers[i], '.html')
    	input_file <- file.path(root_dir, 'OMPARE.Rmd')
    	rmarkdown::render(input = input_file,
    		params = list(topDir = topDir,
    			fusion_method = 'arriba',
        		set_title = set_title,
        		snv_caller = callers[i],
        		tmb = 77.46), 
        	output_dir = output_dir, 
			intermediates_dir = output_dir,
			output_file = output_file)
	}


After running the reports, the project folder will have all output files with plots and tables under ``output`` and all html reports under ``reports``:

.. code-block:: bash

	results/PNOC008-29
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
	│   ├── PNOC008-29_summary_DE_Genes_Down.txt
	│   ├── PNOC008-29_summary_DE_Genes_Up.txt
	│   ├── PNOC008-29_summary_Pathways_Down.txt
	│   ├── PNOC008-29_summary_Pathways_Up.txt
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
	│   ├── PNOC008-29_all.html
	│   ├── PNOC008-29_consensus.html
	│   ├── PNOC008-29_lancet.html
	│   ├── PNOC008-29_mutect2.html
	│   ├── PNOC008-29_strelka2.html
	│   └── PNOC008-29_vardict.html
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
	-p PATIENT, --patient=PATIENT
		Patient Number (1, 2...)

	-w WORKDIR, --workdir=WORKDIR
		OMPARE working directory

	# Example run for PNOC008-21
	Rscript upload_reports.R \
	--patient 21 \
	--wordir /path/to/Projects/OMPARE

Dependencies on HGG-DMG files
=============================

List of files to download and update for each patient:

1. pbta-tgen-gtex-gene-expression-rsem-tpm-collapsed.combined.rds
2. pbta-histologies-base-adapt.tsv
3. pbta-histologies.tsv

Source Location: ``s3://d3b-bix-dev-data-bucket/hgg-dmg-integration/current-data/``
Destination: ``data/reference/pbta``


Dependencies on specific hgg-dmg versions
=========================================

These hgg-dmg files are version dependent:

.. code-block:: bash

	hgg-dmg-integration
	├── 20201109-data
	│   ├── CCDS.bed
	│   ├── StrexomeLite_hg38_liftover_100bp_padded.bed
	│   ├── Strexome_targets_intersect_sorted_padded100.GRCh38.withCCDS.bed
	│   ├── target_cfg.targetcombos.txt
	│   └── xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.withCCDS.bed
	└── 20201202-data
	    ├── CC_based_heatmap_Distance_euclidean_finalLinkage_average_clusterAlg_KM_expct_counts_VST_cluster_and_annotation.tsv
	    ├── pbta-hgat-dx-prog-pm-gene-counts-rsem-expected_count-uncorrected.rds
	    └── pbta-histologies.tsv
