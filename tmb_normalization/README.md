 
## TMB Calculations for TCGA and OpenPBTA Samples Used in OMPARE

Currently, TMB in OMPARE is calculated by dividing the total number of mutations pass filter by the length of a default bed file: `xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed`. This is not the most accurate way since not all samples were processed by the same sequencing panel and hence, would have different effective region of interest (designated by the bed file). This analysis is meant to address this issue.

### Usage
The module does not have a run-script.sh but each script can be run individually as followed:

```
Rscript 00-tcga-not-in-pbta-query-manifest.R
python3 01-download-tcga-mutect2.py
bash 02-prepare-tcga-bed-files.sh
Rscript 03-tcga-sample-kit-bed-anno.R
Rscript 04-calculate-tcga-tmb.R
Rscript 05-calculate-pbta-tmb.R
Rscript 06-tmb-plot.R
```

###### Contents

`00-tcga-not-in-pbta-query-manifest.R` first reads in the current file that we use in OMPARE for TCGA TMB scores (`../references/TCGA_diseasetypes_and_samples_TMBscores.txt`) to identify TCGA projects that are currently in use in OMPARE. The script then uses the project names to query MAF and BAM file manifests from the GDC portal. The results were intersected to get samples that have both MAF and BAM available. Output were saved as `../results/tcga_not_in_pbta_bam_manifest.tsv` and `./results/tcga_not_in_pbta_maf_manifest.tsv`. Additionally, all unique bed files were saved as `../results/tcga_not_in_pbta_unique_bed_file.tsv` and this will be used to query/curl bed files in the following steps.

### input
`../references/TCGA_diseasetypes_and_samples_TMBscores.txt`
### output
`../results/tcga_not_in_pbta_bam_manifest.tsv`
`./results/tcga_not_in_pbta_maf_manifest.tsv`
`../results/tcga_not_in_pbta_unique_bed_file.tsv`


`01-download-tcga-mutect2.py` reads in `../results/tcga_not_in_pbta_maf_manifest.tsv` generated from the previous step and query the GDC API data endpoint for all mutect2 files. Each project (cancer type) has its individual folder. Currently, we need to manually move the mutect2 files from their respective folder to `../../scratch/maf_files` but this can be changed if needed. 

### input
`../results/tcga_not_in_pbta_maf_manifest.tsv`
### output
`../../scratch/maf_files`

`02-prepare-bed-files.sh` reads in `../results/tcga_not_in_pbta_bam_manifest.tsv` that was generated by `00-tcga-not-in-pbta-query-manifest.R`, and query the GDC portal for bed files. **NOTE**: not all files are available for download through the url offered by GDC portal and hence, for the ones that are not available through download, they are downloaded directly from the vendor's website or from this https://mcg.ustc.edu.cn/bsc/cnv/upload/ when not available from the vendor (e.g., earlier versions of the bed file). Additionally, sometimes, TCGA is not sure which bedfile (panel) was used for sequencing - for those cases, the available bed file with the hightest coverage was chosen. The final results were saved as `tbga_not_in_pbta_bed_selected.tsv`. 

`02-prepare-bed-files.sh` also reads in `../results/tcga_in_pbta_bam_manifest.tsv`, which was generated through this previous code available here: https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/tcga-capture-kit-investigation/scripts/get-tcga-capture_kit.py, and query the GDC portal for bed files. 

`02-prepare-bed-files.sh` then used `CrossMap.py` tool to liftOver all bed files from hg18 to hg19 (Gh38). Those bed files include bed files in `tcga_in_pbta`, `tcga_not_in_pbta`, `pbta` and `pnoc008`. 
**NOTE**: bedfiles in `pbta` and `pnoc008` were provided and were not generated by any code in here. The bed files are then sorted and merged and moved to respective folder in the `results` folder for next steps

### input
`../results/tcga_not_in_pbta_bam_manifest.tsv`
`../results/tcga_in_pbta_bam_manifest.tsv`
`../../scratch/pnoc008/ashion_confidential_exome_v2_2nt_pad.bed`

### output
`../results/bed_files/tcga_in_pbta/*.Gh38.bed`
`../results/bed_files/tcga_in_pbta/*.Gh38.bed`
`../references/ashion_confidential_exome_v2_2nt_pad.Gh38.bed`

`03-tcga-sample-kit-bed-anno.R` reads in the sorted, merged and Gh38-corrected bed files from `02-prepare-bed-files.sh` and bam manifest files from `00-tcga-not-in-pbta-query-manifest.R`, calculates the bed length for all bed files and annotates back to the bam manifest with the bed file used and bedlength. 

### input
`../results/bed_files/tcga_in_pbta/*.Gh38.bed`
`../results/bed_files/tcga_in_pbta/*.Gh38.bed`
`../references/ashion_confidential_exome_v2_2nt_pad.Gh38.bed`

### output
`../results/bed_files/tcga_in_pbta/*.Gh38.bed`
`../results/bed_files/tcga_in_pbta/*.Gh38.bed`

`04-calculate-tcga-tmb.R` reads in the mutect2 files for all TCGA projects, the bedfile and bedlength annotated maf files, the unique bed file list with lenghth, and each Gh38 bed files do the TMB calculation.

### input
`../../scratch/maf_files`
`../results/bed_files/{tcga_not_in_pbta,tcga_in_pbta}`
`../results/tcga{_not,}_in_pbta_bam_manifest_with_length.tsv`
`../results/tcga{_not,}_in_pbta_unique_bed_with_length.tsv`

### output
`../results/TCGA_not_in_pbta_diseasetypes_and_samples_TMBscores.new.txt`
`../results/TCGA_in_pbta_diseasetypes_and_samples_TMBscores.new.txt`


`05-calculate-pbta-tmb.R` reads in the histology file of OpenPBTA, the mutect2 file and the default bed files to be used for TMB calculation.

### input
`../references/pbta-histologies.tsv`
`../references/pbta-snv-mutect2.vep.maf.gz`
`../references/pbta_bed_files/xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed`
### output
`../results/PBTA-TMBscores_withdiseastype.txt`


`06-tmb-plot.R` plots individual PNOC008 sample against the TCGA samples in OpenPBTA, TCGA samples not analyzed in OpenPBTA (but previously included in OMPARE) and OpenPBTA samples. At this point, the script is only generating figures for PNO008-32, PNO008-33, and PNO008-34.

### input
`../results/PBTA-TMBscores_withdiseastype.txt`
`../results/TCGA_not_in_pbta_diseasetypes_and_samples_TMBscores.txt`
`../results/TCGA_in_pbta_diseasetypes_and_samples_TMBscores.txt`
`../references/ashion_confidential_exome_v2_2nt_pad.Gh38.bed`

`../references/pnoc008_32.mutect2_somatic.norm.annot.protected.maf`
`../references/pnoc008_33.mutect2_somatic.norm.annot.protected.maf`
`../references/pnoc008_34.mutect2_somatic.norm.annot.protected.maf`

### output
`../plots/PNOC008-32-TMB-plot.png`
`../plots/PNOC008-33-TMB-plot.png`
`../plots/PNOC008-34-TMB-plot.png`


