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

1. R Packages

.. code-block:: bash

	# install packages
	Rscript code/install_pkgs.R

2. Reference Data
   
.. code-block:: bash

	# get reference data from s3
	aws s3 --profile saml s3://d3b-bix-dev-data-bucket/PNOC008/Reference /path/to/OMPARE/data/Reference

Project Organization
====================

1. Clone the OMPARE repository.

2. Download required files from data delivery project:

* Copy Number: `*.CNVs.p.value.txt`
* Copy Number: `*.controlfreec.ratio.txt`
* Expression: `*.genes.results`
* Fusions: `*.arriba.fusions.tsv`
* Fusions: `*.star-fusion.fusion_candidates.final`
* Somatic Variants: `*.maf`
* Germline Variants: `*.hg38_multianno.txt.gz`

3. Organize patient data: 
Run `create_project.R` script to create and organize project folder under data/. This script will also create intermediate folders like `ImmuneScores` and `GSVA` as well as output folders like `Reports` for \*.html reports and `Summary` for excel summary.
   
.. code-block:: bash

	Rscript code/create_project_dir.R --help

	Options:
		-s SOURCEDIR, --sourcedir=SOURCEDIR
			Source directory containing all files from data delivery project

		-d DESTDIR, --destdir=DESTDIR
			Destination directory. Should be /path/to/OMPARE/data/PNOC008-21/ for Patient 13

		-h, --help
			Show this help message and exit

	# Example for Patient PNOC008-21
	Rscript code/create_project.R \
	--sourcedir /path/to/source/PNOC008-21-cavatica-files \
	--destdir /path/to/OMPARE/data/PNOC008-21/

4. Create clinical file using the *create_clinfile.R* script.

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
	--sheet data/Reference/Manifest/PNOC008_Manifest.xlsx \
	--patient PNOC008-21 \
	--dir /path/to/OMPARE/data/PNOC008-21

Steps (3) and (4) should create a folder structure with corresponding files as shown below:

.. code-block:: bash

	# Example for PNOC008-21
	tree /path/to/OMPARE/data/PNOC008-21/
	.
	├── CNV
	│   ├── uuid.controlfreec.CNVs.p.value.txt
	│   └── uuid.controlfreec.ratio.txt
	├── Clinical
	│   └── patient_report.txt
	├── ExpressionGene
	│   └── uuid.rsem.genes.results.gz
	├── Fusions
	│   ├── uuid.STAR.fusion_predictions.abridged.coding_effect.tsv
	│   └── uuid.arriba.fusions.tsv
	├── GSVA
	├── ImmuneScores
	├── MutationsMAF
	│   ├── uuid.consensus_somatic.vep.maf
	│   ├── uuid.gatk.hardfiltered.PASS.vcf.gz.hg38_multianno.txt.gz
	│   ├── uuid.lancet_somatic.vep.maf
	│   ├── uuid.mutect2_somatic.vep.maf
	│   ├── uuid.strelka2_somatic.vep.maf
	│   ├── uuid.vardict_somatic.vep.maf
	├── Reports
	├── Summary

5. Update PNOC008 patient matrices with each new patient data.
   
.. code-block:: bash

	Rscript code/pnoc_format.R

	# Running the script will update the following files:
	data/Reference/PNOC008
	├── PNOC008_TMBscores.rds
	├── PNOC008_TPM_matrix.RDS
	├── PNOC008_clinData.RDS
	├── PNOC008_cnvData_filtered.rds
	├── PNOC008_consensus_mutData_filtered.rds
	├── PNOC008_deg_GTExBrain.rds
	└── PNOC008_fusData_filtered.rds

6. Update GSEA enrichment output with each new patient data.
   
.. code-block:: bash

	Rscript code/gsea_enrichment.R

	# Running the script will update the following files:
	data/Reference/GSEA
	├── PBTA_vs_GTExBrain.RDS
	├── PBTA_vs_PBTA.RDS
	├── PBTA_vs_PBTAHGG.RDS
	├── PNOC008_vs_GTExBrain.RDS
	├── PNOC008_vs_PBTA.RDS
	├── PNOC008_vs_PBTA_HGG.RDS
	├── PNOC008_vs_TCGA_GBM.RDS
	├── TCGA_GBM_vs_GTExBrain.RDS
	└── TCGA_GBM_vs_TCGA_GBM.RDS

7. Excel summary containing up/down pathways and genes of patient of interest vs GTEx Brain, PBTA HGG and PBTA all histologies:

.. code-block:: bash

	Rscript code/tabulate_excel.R --help

	Options:
	-i INPUT, --input=INPUT
		Directory e.g. data/PNOC008-21

	-o OUTPUT, --output=OUTPUT
		output excel file with extension i.e. PNOC008-21_summary.xlsx

	# Example for Patient PNOC008-21
	Rscript code/tabulate_excel.R \
	--input /path/to/OMPARE/data/PNOC008-21 \
	--output PNOC008-21_summary.xlsx

8. Generate markdown report:

.. code-block:: bash

	# topDir is the project directory of current patient
	# fusion_method is the fusion method. Allowed values: *star*, *arriba*, *both* or not specified. (Optional) 
	# set_title is the title for the report. (Optional)
	# snv_pattern is one of the six values for simple variants: *lancet*, *mutect2*, *strelka2*, *vardict*, *consensus*, *all* (all four callers together)
	# tmb (Tumor mutational burden) is set to 77.46.
	setwd(/path/to/OMPARE)
	callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
	for(i in 1:length(callers)) {
	  outputfile <- paste0("data/PNOC008-21/Reports/PNOC008-21_", callers[i], ".html")
	  rmarkdown::render(input = 'OMPARE.Rmd', 
	                    params = list(topDir = 'data/PNOC008-21/',
	                                  fusion_method = 'arriba',
	                                  set_title = 'PNOC008-21 Patient Report',
	                                  snv_pattern = callers[i],
	                                  tmb = 77.46),
	                    output_file = outputfile)
	}


After running step 8, the project folder should have some intermediate and output files:

.. code-block:: bash

	data/PNOC008-21
	├── CNV
	│   ├── uuid.controlfreec.CNVs.p.value.txt
	│   └── uuid.controlfreec.ratio.txt
	├── Clinical
	│   └── patient_report.txt
	├── ExpressionGene
	│   └── uuid.rsem.genes.results.gz
	├── Fusions
	│   ├── uuid.STAR.fusion_predictions.abridged.coding_effect.tsv
	│   └── uuid.arriba.fusions.tsv
	├── GSVA
	│   └── ssgsea_rawScores.txt
	├── ImmuneScores
	│   ├── rawScores_adult.txt
	│   ├── rawScores_pediatric.txt
	│   ├── tisScores.txt
	│   └── topCor_rawScores.txt
	├── MutationsMAF
	│   ├── uuid.lancet_somatic.vep.maf
	│   ├── uuid.mutect2_somatic.vep.maf
	│   ├── uuid.strelka2_somatic.vep.maf
	│   ├── uuid.vardict_somatic.vep.maf
	│   ├── uuid.consensus_somatic.vep.maf
	│   ├── uuid.gatk.PASS.vcf.gz.hg38_multianno.txt.gz
	│   └── mpfDataFormat.txt
	├── Reports
	│   ├── PNOC008-21_all.html
	│   ├── PNOC008-21_consensus.html
	│   ├── PNOC008-21_lancet.html
	│   ├── PNOC008-21_mutect2.html
	│   ├── PNOC008-21_strelka2.html
	│   └── PNOC008-21_vardict.html
	├── Summary
	│   ├── PNOC008-21_summary.xlsx
	│   ├── adultsig_pathways_gen_similar.txt
	│   ├── pbta_pnoc008_umap_output.rds
	│   ├── pediatriccnv_pathways.txt
	│   ├── pediatricsig_pathways_gen_similar.txt
	│   └── tcga_pnoc008_umap_output.rds
	├── complexHeatmap_cgs.png
	├── complexHeatmap_oncogrid.png
	├── complexHeatmap_phgg.png
	└── tmpRCircos.png

Run everything
==============

This single script will take the raw data as input and create output files by:

1. Creating project directory and organize files
2. Creating clinical file
3. Updating PNOC008 data matrices (cnv, mutations, fusions, expression) with each new patient
4. Updating GSEA enrichment outputs with each new patient
5. Generating excel summary
6. Running html reports

.. code-block:: bash
	
	Rscript run_OMPARE.R --help

	Options:
	-p PATIENT, --patient=PATIENT
		Patient Number (1, 2...)

	-s SOURCEDIR, --sourcedir=SOURCEDIR
		Source directory with all files

	-c CLIN_FILE, --clin_file=CLIN_FILE
		PNOC008 Manifest file (.xlsx)

	-w WORKDIR, --workdir=WORKDIR
		OMPARE working directory

	# Example run for PNOC008-21
	Rscript run_OMPARE.R \
	--patient 21 \
	--clin_file data/Reference/Manifest/PNOC008_Manifest.xlsx \
	--workdir ~/Projects/OMPARE \
	--sourcedir ~/Downloads/p21

Upload to data-delivery project
===============================

This script uploads the `Summary/*._summary.xlsx`, `Summary/*._umap_output.rds`, `Reports/*.html` output to the data delivery project folder on cavatica. 

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
	--wordir ~/Projects/OMPARE

