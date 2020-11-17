# Author: Bo Zhang
# Date: 08/10/2020
# Function: Script to generate Oncogrid plot

suppressPackageStartupMessages(library(ComplexHeatmap))

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# oncogrid directory
oncogrid.path <- file.path(ref_dir, 'oncogrid')
oncogrid.path.input <- file.path(oncogrid.path, 'input')
oncogrid.path.output <- file.path(oncogrid.path, 'output')

# matrix
mat = read.table(file.path(oncogrid.path.output, "oncoprint.cohort3.pnoc008.txt"),  header = TRUE, stringsAsFactors=FALSE, sep = "\t",check.names = FALSE)
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat = t(as.matrix(mat))

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#ffffff",col= "#595959"))
  },
  GAI = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#ff4d4d", col = NA))
  },
  LOS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#0D47A1", col = NA))
  },
  MIS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#77b300", col = NA))
  },
  NOS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#80bfff", col = NA))
  },
  FSD = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#1a53ff", col = NA))
  },
  FSI = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#8D6E63", col = NA))
  },
  NOT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#9966ff", col = NA))
  },
  SPS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#E69F00", col = NA))
  },
  IFD = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "#827717", col = NA))
  },
  OVE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.66, gp = gpar(fill = "#dbc6eb", col = NA))
  },
  UNE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.66, gp = gpar(fill = "#709fb0", col = NA))
  },
  FUS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.66, gp = gpar(fill = "#AB47BC", col = NA))
  }
)
col = c("GAI" = "#ff4d4d", "LOS" = "#0D47A1" , "MIS" = "#77b300", "FUS" = "#AB47BC","NOS" ="#80bfff", "FSD" = "#1a53ff", "FSI" ="#8D6E63", "NOT" = "#9966ff","SPS" = "#E69F00","IFD" = "#827717","OVE" = "#dbc6eb","UNE" = "#709fb0")

# read annotation and TMB info
annot_info <- read.table(file.path(oncogrid.path.output, "annotation.cohort3.pnoc008.txt"), header = TRUE, check.names = TRUE)
annot_info <- annot_info %>%
  column_to_rownames('Sample') %>%
  as.data.frame()
tmb_info <- read.table(file.path(oncogrid.path.output, "tmb.cohort3.pnoc008.txt"), header = TRUE, check.names = TRUE)
tmb_info <- tmb_info %>%
  column_to_rownames('Sample') %>%
  as.data.frame()

# annotation 1
ha = HeatmapAnnotation(df = annot_info , col = list(
  Sequencing_Experiment = c("WGS,WXS,RNA-Seq" = "#26734d",
                            "WGS,Panel,RNA-Seq" = "#1a3300",
                            "WGS" = "#009999","RNA-Seq" = "#bac964",
                            "WGS,RNA-Seq"="#5a3d55",
                            "WXS,RNA-Seq"="#96bb7c"),
  Cohort = c("CBTTC" = "#709fb0", 
             "PNOC003" = "#a0c1b8", 
             "PNOC008" = "#f4ebc1",
             "CBTN" = "#a14f62"),
  Tumor_Descriptor = c("Primary" = "#d8b9c3",
                       "Progressive" = "#827397", 
                       "Initial_CNS_Tumor" = "#cee397",
                       "Progressive_Disease_Post_Mortem" = "#4d4c7d",
                       "Recurrence" = "#363062",
                       "Second_Malignancy" = "#005082"),
  Integrated_Diagnosis = c("Brainstem_glioma_Diffuse_intrinsic_pontine_glioma" = "#e1f4f3",
                           "Diffuse_midline_glioma"="#706c61","Gliomatosis_Cerebri" = "#6e5773",
                           "High_grade_glioma" = "#347474"),
  OS_Status = c("DECEASED" = "#e9e2d0",
                "LIVING" = "#7fcd91",
                "NA" = "#f0efef" )),
  gp = gpar(col = "#595959"), simple_anno_size = unit(4, "mm"), annotation_name_side = "left",
  annotation_legend_param = list(
    Sequencing_Experiment = list(nrow =3),
    Cohort = list(nrow=3)
  ))

# split samples
split.samples =read.table(file.path(oncogrid.path.output, "split_samples.cohort3.pnoc008.txt"), header = TRUE, check.names = TRUE)
split.samples = data.frame(split.samples)
subgroup = split.samples$split_samples

# annotation 2
ha1 = HeatmapAnnotation(TMB = anno_barplot(tmb_info$TMB, axis_param = list(side = "right")),
                        annotation_name_side = "left")
amp = ifelse(apply(mat, 1, function(x) sum(grepl("GAI", x) + grepl("LOS",x) + grepl("OVE",x) + grepl("UNE",x))/length(x) > 0),"Copy number alteration and RNAseq Expression","Genetic and Fusion alteration")
amp = factor(amp, levels = c("Genetic and Fusion alteration", "Copy number alteration and RNAseq Expression"))

# gene order
gene_order = scan(file.path(oncogrid.path.input, "sample_cohort3.gene_order"), what = "character")

# oncoprint
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               bottom_annotation = ha1,
               row_names_side = "left",
               pct_side = "right",
               row_split = amp,
               column_split = subgroup,
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD","OVE","UNE"),
                                           labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del","Over Experssion","Under Expression")
               ))


png(filename = file.path(topDir, "output", "complexheatmap_oncogrid.png"), units = "in", width = 25, height = 14, res = 200) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
