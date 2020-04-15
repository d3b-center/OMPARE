# Author: Komal S. Rathi
# Date: 04/13/2020
# Function: Generate patient report

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-p", "--patient"), type = "character",
              help = "Patient Number (1, 2...)"),
  make_option(c("-c", "--clin_file"), type = "character",
              help = "Google sheet link (PNOC008 patients)"),
  make_option(c("-w", "--workdir"), type = "character",
              help = "OMPARE working directory")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
clinical_sheet <- opt$clin_file
p <- opt$patient
workdir <- opt$workdir

# set variables
setwd(workdir)
patient <- paste0('PNOC008-', p)
topDir <- file.path('data', patient)
set_title <- paste0(patient,' Patient Report')
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")

# 1. Create Project directory
cmd1 <- paste0('Rscript create_project_dir.R ', topDir, '/')
system(cmd1)

# 2. Create clinical file
cmd2 <- paste0('Rscript create_clinfile.R -s ', clinical_sheet, ' -p ', patient, ' -d ', topDir)
system(cmd2)

# 3. Update PNOC008 expression matrix for each new patient
cmd3 <- 'Rscript code/pnoc_format.R'
system(cmd3) 

# 4. Run html reports
# fusion_method can be either arriba, star, both or not specified
if(dir.exists(topDir)){
  for(c in 1:length(callers)) {
    outputfile <- paste0(patient,'_',callers[c],'.html')
    outputfile <- file.path(topDir,'Reports',outputfile)
    rmarkdown::render(input = 'OMPARE.Rmd',
                      params = list(topDir = topDir,
                                    fusion_method = 'arriba',
                                    set_title = set_title,
                                    snv_pattern = callers[i],
                                    tmb = 77.46),
                      output_file = outputfile)
  }
} 

# 5. Generate excel summary
cmd5 <- paste0('Rscript code/tabulate_excel.R -i ', topDir, ' -o ', paste0(patient, '_summary.xlsx'))
system(cmd5)
