# Step2: generate patient report

# parameter info:
# fusion_method can be either arriba, star, both or not specified
# fusion_method and set_title is completely optional 
# topDir is mandatory
# snv_pattern can be either of the 6 options: lancet, mutect2, strelka2, vardict, consensus or all

# PNOC008-04
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
for(i in 1:length(callers)) {
  outputfile <- paste0("data/PNOC008-04/Reports/PNOC008_04_", callers[i], ".html")
  rmarkdown::render(input = 'OMPARE.Rmd', 
                    params = list(topDir = 'data/PNOC008-04/',
                                  fusion_method = 'arriba',
                                  set_title = 'PNOC008-04 Patient Report',
                                  snv_pattern = callers[i],
                                  tmb = 77.46),
                    output_file = outputfile)
}

# PNOC008-02
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
for(i in 1:length(callers)) {
  outputfile <- paste0("data/PNOC008-02/Reports/PNOC008_02_", callers[i], ".html")
  rmarkdown::render(input = 'OMPARE.Rmd', 
                    params = list(topDir = 'data/PNOC008-02/',
                                  fusion_method = 'arriba',
                                  set_title = 'PNOC008-02 Patient Report',
                                  snv_pattern = callers[i],
                                  tmb = 77.46),
                    output_file = outputfile)
}
