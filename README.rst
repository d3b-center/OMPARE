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
	# currently data is available in the repo. In the future, we will pull from s3
	aws s3 --profile saml s3://d3b-bix-dev-data-bucket/PNOC008/Reference /path/to/OMPARE/data/Reference

Project Organization
====================

1. Clone the OMPARE repository.

2. Download required files from data delivery project:

* Copy Number: .CNVs.p.value.txt
* Copy Number: .controlfreec.ratio.txt
* Expression: .genes.results
* Fusions: .arriba.fusions.tsv
* Fusions: .star-fusion.fusion_candidates.final
* Somatic Variants: .maf
* Germline Variants: .hg38_multianno.txt.gz

3. Organize patient data: 
Run *create_project.R* script to create and organize project folder under data/. This script will also create intermediate folders like *ImmuneScores* and *GSVA* as well as output folders like *Reports* for .html reports and *Summary* for excel summary.
   
.. code-block:: bash

	# Run script to create project folder
	# -s source directory will data dump from cavatica
	# -d destination directory. Should be /path/to/OMPARE/data/subjectID
	# Example for Patient PNOC008-13
	Rscript create_project.R -s /path/to/source/PNOC008-13-cavatica-files -d /path/to/OMPARE/data/PNOC008-13/

4. Create clinical file using the *create_clinfile.R* script.

.. code-block:: bash

	# Run script to create clinical file
	# -s is the env variable PNOC008_MANIFEST which is the link to the manifest on google sheets
	# -p parameter should match the subjectID in the manifest so check that before running
	# -d is the path to project directory. Should be /path/to/OMPARE/data/subjectID
	# Example for Patient PNOC008-13
	Rscript create_clinfile.R -s $PNOC008_MANIFEST -p PNOC008-8 -d /path/to/OMPARE/data/PNOC008-13

Steps (3) and (4) should create a folder structure as shown below:

.. code-block:: bash

	# Example for PNOC008-13
	tree /path/to/OMPARE/data/PNOC008-13/
	.
	├── CNV
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.controlfreec.CNVs.p.value.txt
	│   └── 86319bd2-424d-44bf-928b-03a0ad97ee54.controlfreec.ratio.txt
	├── Clinical
	│   └── patient_report.txt
	├── ExpressionGene
	│   └── d27a8113-7acb-4303-a018-e1383c2673af.rsem.genes.results.gz
	├── Fusions
	│   ├── d27a8113-7acb-4303-a018-e1383c2673af.STAR.fusion_predictions.abridged.coding_effect.tsv
	│   └── d27a8113-7acb-4303-a018-e1383c2673af.arriba.fusions.tsv
	├── GSVA
	├── ImmuneScores
	├── MutationsMAF
	│   ├── 5ff3b001-bbdd-4640-829e-5d340aa90482.consensus_somatic.vep.maf
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.gatk.hardfiltered.PASS.vcf.gz.hg38_multianno.txt.gz
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.lancet_somatic.vep.maf
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.mutect2_somatic.vep.maf
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.strelka2_somatic.vep.maf
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.vardict_somatic.vep.maf
	├── Reports
	├── Summary

5. Update PNOC008 expression matrix for each new patient. This will update PNOC008_matrix.RDS (expression matrix) and PNOC008_clinData.RDS (clinical file) under /path/to/OMPARE/data/Reference/PNOC008.
   
.. code-block:: bash

	Rscript code/pnoc_format.R   

6. Report:

.. code-block:: bash

	# topDir is your project directory
	# fusion_method is the fusion method. Allowed values: *star*, *arriba*, *both* or not specified. (Optional) 
	# set_title is the title for the report. (Optional)
	# snv_pattern is one of the six values for simple variants: *lancet*, *mutect2*, *strelka2*, *vardict*, *consensus*, *all* (all four callers together)
	# tmb (Tumor mutational burden) is set to 77.46.
	setwd(/path/to/OMPARE)
	callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
	for(i in 1:length(callers)) {
	  outputfile <- paste0("data/PNOC008-13/Reports/PNOC008-13_", callers[i], ".html")
	  rmarkdown::render(input = 'OMPARE.Rmd', 
	                    params = list(topDir = 'data/PNOC008-13/',
	                                  fusion_method = 'arriba',
	                                  set_title = 'PNOC008-13 Patient Report',
	                                  snv_pattern = callers[i],
	                                  tmb = 77.46),
	                    output_file = outputfile)
	}

7. Excel summary:

.. code-block:: bash

	# create excel summary file
	Rscript code/tabulate_excel.R -i /path/to/OMPARE/data/PNOC008-13 -o PNOC008-13_summary.xlsx

After running step 7, the project folder should have some intermediate and output files:

* tmpRCircos.png: circos plot 
* ImmuneScores/rawScores.txt: raw scores from immune profile function
* GSVA/hallmark_rawScores.txt: ssGSEA raw scores
* Reports/\*.html for each individual caller, consensus and all callers together.
* Summary/\*.excel summary report.

.. code-block:: bash

	tree data/PNOC008-13/
	.
	├── CNV
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.controlfreec.CNVs.p.value.txt
	│   └── 86319bd2-424d-44bf-928b-03a0ad97ee54.controlfreec.ratio.txt
	├── Clinical
	│   └── patient_report.txt
	├── ExpressionGene
	│   └── d27a8113-7acb-4303-a018-e1383c2673af.rsem.genes.results.gz
	├── Fusions
	│   ├── d27a8113-7acb-4303-a018-e1383c2673af.STAR.fusion_predictions.abridged.coding_effect.tsv
	│   └── d27a8113-7acb-4303-a018-e1383c2673af.arriba.fusions.tsv
	├── GSVA
	│   └── hallmark_rawScores.txt
	├── ImmuneScores
	│   ├── rawScores.txt
	│   └── topCor_rawScores.txt
	├── MutationsMAF
	│   ├── 5ff3b001-bbdd-4640-829e-5d340aa90482.consensus_somatic.vep.maf
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.gatk.hardfiltered.PASS.vcf.gz.hg38_multianno.txt.gz
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.lancet_somatic.vep.maf
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.mutect2_somatic.vep.maf
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.strelka2_somatic.vep.maf
	│   ├── 86319bd2-424d-44bf-928b-03a0ad97ee54.vardict_somatic.vep.maf
	│   └── mpfDataFormat.txt
	├── Reports
	│   ├── PNOC008_13_all.html
	│   ├── PNOC008_13_consensus.html
	│   ├── PNOC008_13_lancet.html
	│   ├── PNOC008_13_mutect2.html
	│   ├── PNOC008_13_strelka2.html
	│   └── PNOC008_13_vardict.html
	├── Summary
	│   ├── PNOC008-13_summary.xlsx
	│   └── up_pathways_gen_similar.txt
	└── tmpRCircos.png


Run everything
==============

This single script will take the raw data as input and create output files by:

1. Creating project directory and organize files
2. Creating clinical file
3. Updating PNOC008 expression matrix for each new patient
4. Running html reports
5. Generating excel summary

.. code-block:: bash
	
	Rscript run_OMPARE.R -p 13 -c <link_to_google_sheet> -w <OMPARE_directory>
