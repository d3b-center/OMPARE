# script to copy output files generated after running scripts 00-05 to data/reference/tmb folder

# set root directory and other directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
tmb_results_dir <- file.path(root_dir, 'tmb_normalization', 'results')
tmb_ref_dir <- file.path(root_dir, 'data', 'reference', 'tmb')

# cp output files generated after running scripts 00-05 to data/reference/tmb folder
pbta <- paste('cp', 
      file.path(tmb_results_dir, 'PBTA-TMBscores_withdiseastype.txt'), 
      file.path(tmb_ref_dir, 'PBTA-TMBscores_withdiseastype.txt'))
system(pbta)

tcga_not_in_pbta <- paste('cp', 
             file.path(tmb_results_dir, 'TCGA_not_in_pbta_diseasetypes_and_samples_TMBscores.txt'), 
             file.path(tmb_ref_dir, 'TCGA_not_in_pbta_diseasetypes_and_samples_TMBscores.txt'))
system(tcga_not_in_pbta)

tcga_in_pbta <- paste('cp', 
             file.path(tmb_results_dir, 'TCGA_in_pbta_diseasetypes_and_samples_TMBscores.txt'), 
             file.path(tmb_ref_dir, 'TCGA_in_pbta_diseasetypes_and_samples_TMBscores.txt'))
system(tcga_in_pbta)