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

3. Download reference Data:
   
.. code-block:: bash

	# get reference data from s3
	aws s3 --profile saml s3://d3b-bix-dev-data-bucket/PNOC008/Reference /path/to/OMPARE/data/reference

4. Download input files from data delivery project:

* Copy Number: ``*.CNVs.p.value.txt``
* Copy Number: ``*.controlfreec.ratio.txt``
* Expression: ``*.genes.results``
* Fusions: ``*.arriba.fusions.tsv``
* Fusions: ``*.star-fusion.fusion_candidates.final``
* Somatic Variants: ``*.maf``
* Germline Variants: ``*.hg38_multianno.txt.gz``

Scripts
=======

Master script
-------------

**run_OMPARE.R**: Master script that runs the following scripts:
   
1. **code/create_project_dir.R**: creating project directory and organize files.
2. **code/create_clinfile.R**: creating clinical file for patint of interest.
3. **code/patient_level_analyses/pnoc_format.R**: updating PNOC008 data matrices (cnv, mutations, fusions, expression) with each new patient.
4. **code/patient_level_analyses/gsea_enrichment.R**: updating GSEA enrichment outputs with each new patient.
5. **code/patient_level_analyses/tabulate_excel.R**: generating excel file with up/down pathways and genes.
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

**code/create_project_dir.**R: this script creates and organizes input files under ``results``. Creates ``output`` folder to store all output for plots and tables reported and ``reports`` folder to store all html output.
   
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

Update PNOC008 data matrices:
-----------------------------

**code/patient_level_analyses/pnoc_format.R**: this script updates the 008 patient matrices (cnv, mutations, fusions, expression) by adding current patient of interest
   
.. code-block:: bash

	Rscript code/patient_level_analyses/pnoc_format.R

	# Running the script will update the following files:
	data/reference/PNOC008
	├── PNOC008_TMBscores.rds
	├── PNOC008_TPM_matrix.RDS
	├── PNOC008_clinData.RDS
	├── PNOC008_cnvData_filtered.rds
	├── PNOC008_consensus_mutData_filtered.rds
	├── PNOC008_deg_GTExBrain.rds
	└── PNOC008_fusData_filtered.rds

Update GSEA enrichment:
-----------------------

**code/patient_level_analyses/gsea_enrichment.R**: this script will update GSEA enrichment output with each new patient data.
   
.. code-block:: bash

	Rscript code/patient_level_analyses/gsea_enrichment.R

	# Running the script will update the following files:
	data/reference/GSEA
	├── PBTA_vs_GTExBrain.RDS
	├── PBTA_vs_PBTA.RDS
	├── PBTA_vs_PBTAHGG.RDS
	├── PNOC008_vs_GTExBrain.RDS
	├── PNOC008_vs_PBTA.RDS
	├── PNOC008_vs_PBTA_HGG.RDS
	├── PNOC008_vs_TCGA_GBM.RDS
	├── TCGA_GBM_vs_GTExBrain.RDS
	└── TCGA_GBM_vs_TCGA_GBM.RDS

Excel file with differential results:
-------------------------------------

**code/patient_level_analyses/tabulate_excel.R**: this script will create an excel summary containing up/down pathways and genes of patient of interest vs ``GTEx Brain``, ``PBTA HGG`` and ``PBTA all histologies``:

.. code-block:: bash

	Rscript code/patient_level_analyses/tabulate_excel.R --help

	Options:
	-i INPUT, --input=INPUT
		Directory e.g. results/PNOC008-21

	-o OUTPUT, --output=OUTPUT
		output excel file with extension i.e. PNOC008-21_summary.xlsx

	# Example for Patient PNOC008-21
	Rscript code/tabulate_excel.R \
	--input /path/to/OMPARE/results/PNOC008-21 \
	--output PNOC008-21_summary.xlsx

HTML reports:
-------------

8. Generate markdown report:

.. code-block:: bash

	# topDir is the project directory of current patient
	# fusion_method is the fusion method. Allowed values: star, arriba, both or not specified. (Optional) 
	# set_title is the title for the report. (Optional)
	# snv_pattern is one of the six values for simple variants: lancet, mutect2, strelka2, vardict, consensus, all (all four callers together)
	# tmb (Tumor mutational burden) is set to 77.46.
	setwd(/path/to/OMPARE)
	callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
	for(i in 1:length(callers)) {
	  outputfile <- paste0("results/PNOC008-21/Reports/PNOC008-21_", callers[i], ".html")
	  rmarkdown::render(input = 'OMPARE.Rmd', 
	                    params = list(topDir = 'results/PNOC008-21/',
	                                  fusion_method = 'arriba',
	                                  set_title = 'PNOC008-21 Patient Report',
	                                  snv_pattern = callers[i],
	                                  tmb = 77.46),
	                    output_file = outputfile)
	}


After running the reports, the project folder will have all output files with plots and tables under ``output`` and all html reports under ``reports``:

.. code-block:: bash

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
	│   ├── PNOC008-21_summary.xlsx
	│   ├── adult_immune_profile.rds
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
	│   ├── filtered_germline_vars.rds
	│   ├── immune_scores_adult.txt
	│   ├── immune_scores_pediatric.txt
	│   ├── immune_scores_topcor_pediatric.txt
	│   ├── kaplan_meier_adult.rds
	│   ├── kaplan_meier_pediatric.rds
	│   ├── mutational_analysis_pediatric.rds
	│   ├── network_plot_output.rds
	│   ├── pathway_analysis_adult.rds
	│   ├── pathway_analysis_pediatric.rds
	│   ├── pbta_pnoc008_umap_output.rds
	│   ├── pediatric_immune_profile.rds
	│   ├── pediatric_topcor_immune_profile.rds
	│   ├── ssgsea_pediatric.rds
	│   ├── ssgsea_scores_pediatric.txt
	│   ├── tcga_pnoc008_umap_output.rds
	│   ├── tis_profile.rds
	│   ├── tis_scores.txt
	│   ├── tmb_profile_output.rds
	│   ├── transciptomically_similar_adult.rds
	│   ├── transciptomically_similar_pediatric.rds
	│   └── tumor_signature_output.rds
	├── reports
	│   ├── PNOC008-21_all.html
	│   ├── PNOC008-21_consensus.html
	│   ├── PNOC008-21_lancet.html
	│   ├── PNOC008-21_mutect2.html
	│   ├── PNOC008-21_strelka2.html
	│   └── PNOC008-21_vardict.html
	└── simple-variants
	    ├── uuid.lancet_somatic.vep.maf
	    ├── uuid.mutect2_somatic.vep.maf
	    ├── uuid.strelka2_somatic.vep.maf
	    ├── uuid.vardict_somatic.vep.maf
	    ├── uuid.consensus_somatic.vep.maf
	    └── uuid.gatk.PASS.vcf.gz.hg38_multianno.txt.gz

Upload to data-delivery project
-------------------------------

**upload_reports.R**: this script uploads the files under ``reports`` and ``output`` folder to the data delivery project folder on cavatica. 

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

