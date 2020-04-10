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

# PNOC008-06
# reports
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
for(i in 1:length(callers)) {
  outputfile <- paste0("data/PNOC008-06/Reports/PNOC008_06_", callers[i], ".html")
  rmarkdown::render(input = 'OMPARE.Rmd', 
                    params = list(topDir = 'data/PNOC008-06/',
                                  fusion_method = 'arriba',
                                  set_title = 'PNOC008-06 Patient Report',
                                  snv_pattern = callers[i],
                                  tmb = 77.46),
                    output_file = outputfile)
}
# summary
system("Rscript code/tabulate_excel.R -i data/PNOC008-06 -o PNOC008-06_summary-v2.xlsx")

# PNOC008-05-NANT
# reports
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
for(i in 1:length(callers)) {
  outputfile <- paste0("data/PNOC008-05-NANT/Reports/PNOC008_05_NANT_", callers[i], ".html")
  rmarkdown::render(input = 'OMPARE.Rmd', 
                    params = list(topDir = 'data/PNOC008-05-NANT/',
                                  fusion_method = 'arriba',
                                  set_title = 'PNOC008-05-NANT Patient Report',
                                  snv_pattern = callers[i],
                                  tmb = 77.46),
                    output_file = outputfile)
}
# summary
system("Rscript code/tabulate_excel.R -i data/PNOC008-05-NANT -o PNOC008-05-NANT_summary.xlsx")

# PNOC008-05-CHOP
# reports
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
for(i in 1:length(callers)) {
  outputfile <- paste0("data/PNOC008-05-CHOP/Reports/PNOC008_05_CHOP_", callers[i], ".html")
  rmarkdown::render(input = 'OMPARE.Rmd', 
                    params = list(topDir = 'data/PNOC008-05-CHOP/',
                                  fusion_method = 'arriba',
                                  set_title = 'PNOC008-05-CHOP Patient Report',
                                  snv_pattern = callers[i],
                                  tmb = 77.46),
                    output_file = outputfile)
}
# summary
system("Rscript code/tabulate_excel.R -i data/PNOC008-05-CHOP -o PNOC008-05-CHOP_summary.xlsx")

# PNOC008-08
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

# PNOC008-09
# reports
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
for(i in 1:length(callers)) {
  outputfile <- paste0("data/PNOC008-09/Reports/PNOC008_09_", callers[i], ".html")
  rmarkdown::render(input = 'OMPARE.Rmd', 
                    params = list(topDir = 'data/PNOC008-09/',
                                  fusion_method = 'arriba',
                                  set_title = 'PNOC008-09 Patient Report',
                                  snv_pattern = callers[i],
                                  tmb = 77.46),
                    output_file = outputfile)
}
# summary
system("Rscript code/tabulate_excel.R -i data/PNOC008-09 -o PNOC008-09_summary.xlsx")

# To run all reports
# update PNOC008 expression matrix for each new patient
system('Rscript code/pnoc_format.R') 
patients <- 15
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")
for(p in 1:patients){
  topDir <- paste0('data/PNOC008-',p,'/')
  set_title <- paste0('PNOC008-',p,' Patient Report')
  if(dir.exists(topDir)){
    for(c in 1:length(callers)) {
      outputfile <- paste0('data/PNOC008-',p,'/Reports/PNOC008_',p,'_', callers[c], '.html')
      rmarkdown::render(input = 'OMPARE.Rmd',
                        params = list(topDir = topDir,
                                      fusion_method = 'arriba',
                                      set_title = set_title,
                                      snv_pattern = callers[i],
                                      tmb = 77.46),
                        output_file = outputfile)
    }
  }
}