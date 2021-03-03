# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'diffreg_pathways_barplot.R'))

# input files
pathways_up <- read.delim(file.path(topDir, "output", paste0(sampleInfo$subjectID, '_summary_Pathways_Up.txt')))
pathways_down <- read.delim(file.path(topDir, "output", paste0(sampleInfo$subjectID, '_summary_Pathways_Down.txt')))

# call function
fname <- file.path(topDir, "output", "diffreg_pathways_barplot_output.rds")
diffreg_pathways_barplot_gtex <- diffreg_pathways_barplot(pathways_up, pathways_down, comparison_study = 'GTExBrain_1152')
diffreg_pathways_barplot_pbta_hgg <- diffreg_pathways_barplot(pathways_up, pathways_down, comparison_study = 'PBTA_HGG_182')
diffreg_pathways_barplot_pbta <- diffreg_pathways_barplot(pathways_up, pathways_down, comparison_study = 'PBTA_All_1028')
diffreg_pathways_barplot_output <- list(diffreg_pathways_barplot_gtex = diffreg_pathways_barplot_gtex,
     diffreg_pathways_barplot_pbta_hgg = diffreg_pathways_barplot_pbta_hgg,
     diffreg_pathways_barplot_pbta = diffreg_pathways_barplot_pbta)
saveRDS(diffreg_pathways_barplot_output, file = fname)
