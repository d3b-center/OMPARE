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

Input files
===========

Instructions:
	
- Create project folder under *data/*. 
- Keep the subdirectory names and file extensions consistent.

.. code-block:: bash

	# Directory structure:
	data/PNOC008
	├── CNV
	│   └── 8a03c927-7f37-4605-ba6e-73a885cade6c.CNVs
	├── Clinical
	│   └── patient_report.txt
	├── ExpressionGene
	│   └── 7316-903_577716.genes.results
	├── Fusions
	│   ├── 7316-535.local.transcript.converted.pe.star-fusion.fusion_candidates.final
	│   └── eb06d52a-110b-4c0e-9e90-94ea68d2f698.arriba.fusions.tsv
	├── ImmuneScores
	│   └── rawScores.txt
	├── MutationsMAF
	│   ├── 5c3eace5-950a-4a05-81ed-5c04b4a0a367.strelka.vep.maf
	│   ├── 8a03c927-7f37-4605-ba6e-73a885cade6c.strelka.vep.maf
	│   └── ae4ce725-d3c4-455d-822e-6c5067444b5e.strelka.vep.maf
	└── tmpRCircos.png

- Files provided by user:

    + CNV/\*.CNVs (Optional)
    + Clinical/patient_report.txt (Optional)
    + ExpressionGene/\*.genes.results (Required)
    + Fusions/\*.arriba.fusions.tsv (Optional)
    + Fusions/\*.star-fusion.fusion_candidates.final (Optional)
    + MutationsMAF/\*.maf (Optional)

- Files created upon execution:

    + *tmpRCircos.png*. Requires Fusion data. 
    + *ImmuneScores/rawScores.txt*. Requires Expression data.

Running the code
================

Input Parameters: 

- *topDir* is your project directory. (Required)
- *fusion_method* is the fusion method, currently only one is used. Allowed values: *star* or *arriba*. (Optional) 
- *set_title* is the title for the report. (Optional)

.. code-block:: bash

	# run with custom parameter values
	rmarkdown::render(input = 'OMPARE.Rmd', 
                  	  params = list(topDir = 'data/PNOC008/', 
                  	  fusion_method = 'arriba',
                  	  set_title = 'PNOC008 Report'))

