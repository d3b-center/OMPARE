# Function: script to generate oncogrid plot

suppressPackageStartupMessages({
  library(ComplexHeatmap)
})

plot_oncogrid <- function(input_dir, output_file, width, height){
  # matrix
  mat = readRDS(file.path(input_dir, "oncogrid_input_matrix.rds"))
  mat[is.na(mat)] = ""
  mat[mat == "NA"] <- ""
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
  annot_info <- readRDS(file.path(input_dir, "oncogrid_input_annotation.rds"))
  annot_info <- annot_info %>%
    mutate(Tumor_Descriptor = gsub(" ", "_", Tumor_Descriptor)) %>%
    remove_rownames() %>%
    column_to_rownames('Sample') %>%
    dplyr::select(-c(cohort_participant_id)) %>%
    as.data.frame()
  annot_info[is.na(annot_info)] <- "NA"
  
  # annotation 1
  ha = HeatmapAnnotation(df = annot_info , col = list(
    Cohort = c("POI" = "#f4ebc1",
               "PBTA" = "#a14f62",
               "PNOC" = "#8656b3",
               "DGD" = "#f1c232",
               "CBTN" = "#2986cc",
               "GMKF" = "#6aa84f",
               "TARGET" = "#f44336"),
    Tumor_Descriptor = c("Primary" = "#d8b9c3",
                         "Tumor" = "#d8b9c3",
                         "Progressive" = "#827397", 
                         "Initial_CNS_Tumor" = "#cee397",
                         "Progressive_Disease_Post_Mortem" = "#4d4c7d",
                         "Recurrence" = "#363062",
                         "Second_Malignancy" = "#005082",
                         "NA" = "#f0efef",
                         "Not_Applicable" = "#f0efef"),
    OS_Status = c("DECEASED" = "#e9e2d0",
                  "LIVING" = "#7fcd91",
                  "NA" = "#f0efef")),
    gp = gpar(col = "#595959"), simple_anno_size = unit(4, "mm"), annotation_name_side = "left",
    annotation_legend_param = list(
      Sequencing_Experiment = list(nrow = 4),
      Cohort = list(nrow = 3)
    ))
  
  
  # annotation 2
  amp = ifelse(apply(mat, 1, function(x) sum(grepl("GAI", x) + grepl("LOS",x) + grepl("OVE",x) + grepl("UNE",x))/length(x) > 0),"Copy number alteration and RNAseq Expression","Genetic and Fusion alteration")
  amp = factor(amp, levels = c("Genetic and Fusion alteration", "Copy number alteration and RNAseq Expression"))
  
  # oncoprint
  ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
                 alter_fun = alter_fun, col = col, show_column_names =TRUE, column_names_gp = gpar(fontsize = 9),
                 column_names_side = "top",
                 top_annotation = ha,
                 right_annotation = NULL,
                 row_names_side = "left",
                 pct_side = "right",
                 row_split = amp,
                 remove_empty_rows = TRUE,
                 heatmap_legend_param = list(title = "Alterations", nrow = 9, title_position = "topleft", direction = "horizontal",
                                             at = c("GAI","LOS","MIS","FUS","NOS","FSD","FSI","SPS","IFD","OVE","UNE"),
                                             labels = c("Copy gain", "Copy loss", "Misense","Gene Fusion","Nonsense","Frame_Shift_Del","Frame_Shift_Ins","Splice site","In_Frame_Del","Over Experssion","Under Expression")
                 ))
  
  pdf(file = output_file, width = width, height = height) 
  draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}
