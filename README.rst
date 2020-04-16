.. |date| date::

********************
Omics Patient Report
********************

:authors: Pichai Raman, Komal S Rathi
:contact: ramanp@email.chop.edu
:organization: D3B, CHOP
:status: This is "work in progress"
:date: |date|

.. meta::
   :keywords: omics, report, flexboard, 2019
   :description: Omics Patient Report

Prerequisites
=============

.. code-block:: bash

	# install packages
	Rscript code/install_pkgs.R

Project Folder Organization
===========================

**Input files:**

* Copy Number: CNV/\*.CNVs.p.value.txt (Optional)
* Clinical: Clinical/patient_report.txt (Optional)
* Expression: ExpressionGene/\*.genes.results (Optional)
* Fusions: Fusions/\*.arriba.fusions.tsv (Optional)
* Fusions: Fusions/\*.star-fusion.fusion_candidates.final (Optional)
* Somatic Variants: MutationsMAF/\*.maf (Optional)
* Germline Variants: MutationsMAF/\*.hg38_multianno.txt.gz (Optional)


**Instructions:**
	
1. Clone this repo.
2. Create project folder using the convention *OMPARE/data/subjectID*.
3. Dump all downloaded data under the project folder.
4. Run *create_project.R* script to create and organize project folder. This script will also create intermediate folders like *ImmuneScores* and output folders like *Reports* for .html reports and *Summary* for excel summary.
5. *patient_report.txt* can either be created manually or using the *create_clinfile.R* script.

.. code-block:: bash

	# Run script to organize files into corresponding folders
	# this just needs one argument: path to project directory

	Rscript create_project.R data/PNOC008-08/

	# Run script to create clinical file
	# -s is the env variable PNOC008_MANIFEST which is the link to the manifest on google sheets
	# -p parameter should match the subjectID in the manifest so check that before running
	# -d is the path to project directory

	Rscript create_clinfile.R -s $PNOC008_MANIFEST -p PNOC008-8 -d data/PNOC008-08

**Output files:**

The above scripts should create a folder structure as shown below:

.. code-block:: bash

	tree data/PNOC008-08/
	.
	├── CNV
	│   └── f12011c0-2981-4d54-9678-79988d67ded8.controlfreec.CNVs.p.value.txt
	├── Clinical
	│   └── patient_report.txt
	├── ExpressionGene
	│   └── bf7cafc7-9f33-4d6a-a088-e1794a731232.rsem.genes.results.gz
	├── Fusions
	│   ├── bf7cafc7-9f33-4d6a-a088-e1794a731232.STAR.fusion_predictions.abridged.coding_effect.tsv
	│   └── bf7cafc7-9f33-4d6a-a088-e1794a731232.arriba.fusions.tsv
	├── ImmuneScores
	├── MutationsMAF
	│   ├── 5681def8-e594-4866-b612-26ad07a8f20b.gatk.hardfiltered.PASS.vcf.gz.hg38_multianno.txt.gz
	│   ├── bf2a2a9f-0251-4017-aaeb-0c3b97690273.consensus_somatic.vep.maf
	│   ├── f12011c0-2981-4d54-9678-79988d67ded8.lancet_somatic.vep.maf
	│   ├── f12011c0-2981-4d54-9678-79988d67ded8.mutect2_somatic.vep.maf
	│   ├── f12011c0-2981-4d54-9678-79988d67ded8.strelka2_somatic.vep.maf
	│   ├── f12011c0-2981-4d54-9678-79988d67ded8.vardict_somatic.vep.maf
	├── Reports
	└── Summary


Report Generation
=================

**Input Parameters:** 

- *topDir* is your project directory. (Required)
- *fusion_method* is the fusion method. Allowed values: *star*, *arriba*, *both* or not specified. (Optional) 
- *set_title* is the title for the report. (Optional)
- *snv_pattern* is one of the six values for simple variants: *lancet*, *mutect2*, *strelka2*, *vardict*, *consensus*, *all* (all four callers together)
- *tmb* (Tumor mutational burden) is set to 77.46.
  
**NOTE**: Easiest way to run the report is the use the template below and just replace the `subjectID`.


**Instructions:**

.. code-block:: bash

	# e.g. of run using PNOC008-08
	setwd('/path/to/OMPARE/')
	# reports
	callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
	for(i in 1:length(callers)) {
	  outputfile <- paste0("data/PNOC008-08/Reports/PNOC008_08_", callers[i], ".html")
	  rmarkdown::render(input = 'OMPARE.Rmd', 
	                    params = list(topDir = 'data/PNOC008-08/',
	                                  fusion_method = 'arriba',
	                                  set_title = 'PNOC008-08 Patient Report',
	                                  snv_pattern = callers[i],
	                                  tmb = 77.46),
	                    output_file = outputfile)
	}
	# summary
	system("Rscript code/tabulate_excel.R -i data/PNOC008-08 -o PNOC008-08_summary.xlsx")


**Output files:**

These are some intermediate and final files created after running the code:

* tmpRCircos.png: Requires Fusion data. 
* ImmuneScores/rawScores.txt: Requires Expression data.
* Reports/\*.html for each individual caller, consensus and all callers together.
* Summary/\*.excel summary report.

The project folder will look like this:

.. code-block:: bash

	tree data/PNOC008-08/
	.
	├── CNV
	│   └── f12011c0-2981-4d54-9678-79988d67ded8.controlfreec.CNVs.p.value.txt
	├── Clinical
	│   └── patient_report.txt
	├── ExpressionGene
	│   └── bf7cafc7-9f33-4d6a-a088-e1794a731232.rsem.genes.results.gz
	├── Fusions
	│   ├── bf7cafc7-9f33-4d6a-a088-e1794a731232.STAR.fusion_predictions.abridged.coding_effect.tsv
	│   └── bf7cafc7-9f33-4d6a-a088-e1794a731232.arriba.fusions.tsv
	├── ImmuneScores
	│   └── rawScores.txt
	├── MutationsMAF
	│   ├── 5681def8-e594-4866-b612-26ad07a8f20b.gatk.hardfiltered.PASS.vcf.gz.hg38_multianno.txt.gz
	│   ├── bf2a2a9f-0251-4017-aaeb-0c3b97690273.consensus_somatic.vep.maf
	│   ├── f12011c0-2981-4d54-9678-79988d67ded8.lancet_somatic.vep.maf
	│   ├── f12011c0-2981-4d54-9678-79988d67ded8.mutect2_somatic.vep.maf
	│   ├── f12011c0-2981-4d54-9678-79988d67ded8.strelka2_somatic.vep.maf
	│   ├── f12011c0-2981-4d54-9678-79988d67ded8.vardict_somatic.vep.maf
	│   └── mpfDataFormat.txt
	├── Reports
	│   ├── PNOC008_08_all.html
	│   ├── PNOC008_08_consensus.html
	│   ├── PNOC008_08_lancet.html
	│   ├── PNOC008_08_mutect2.html
	│   ├── PNOC008_08_strelka2.html
	│   └── PNOC008_08_vardict.html
	└── Summary
	    └── PNOC008-08_summary.xlsx


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
