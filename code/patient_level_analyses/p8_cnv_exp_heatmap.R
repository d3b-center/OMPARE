# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'create_heatmap.R'))

# heatmap with phgg genes
create_heatmap(fname = file.path(topDir, "output", "complexheatmap_phgg.png"),
               genelist = genelist.heatmap$pHGG_Gene_List,
               plot.layout = "h")

# heatmap with cgs genes
create_heatmap(fname = file.path(topDir, "output", "complexheatmap_cgs.png"),
               genelist = genelist.heatmap$Cancer_Gene_Census_CNVs,
               plot.layout = "h")
