# Step1: code to create and organize directory structure
# Rscript create_project <path to project directory with all files>

# required
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("At most one argument must be supplied", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  print("Proceed")
}

topDir <- args[1]

# specify directories
cnvdir <- paste0(topDir, 'CNV')
clinicaldir <- paste0(topDir, 'Clinical')
exprdir <- paste0(topDir, 'ExpressionGene')
fusionsdir <- paste0(topDir, 'Fusions')
immunescores <- paste0(topDir, 'ImmuneScores')
mutdir <- paste0(topDir, 'MutationsMAF')
reports <- paste0(topDir, 'Reports')
summary <- paste0(topDir, 'Summary')

# create directories
system(paste0('mkdir -p ', cnvdir))
system(paste0('mkdir -p ', clinicaldir))
system(paste0('mkdir -p ', exprdir))
system(paste0('mkdir -p ', fusionsdir))
system(paste0('mkdir -p ', immunescores))
system(paste0('mkdir -p ', mutdir))
system(paste0('mkdir -p ', reports))
system(paste0('mkdir -p ', summary))

# organize data
cmd <- paste0('mv ', topDir, '*CNVs* ', cnvdir)
system(cmd)
cmd <- paste0('mv ', topDir, 'patient_report.txt ', clinicaldir)
system(cmd)
cmd <- paste0('mv ', topDir, '*rsem* ', exprdir)
system(cmd)
cmd <- paste0('mv ', topDir, '*fusion* ', fusionsdir)
system(cmd)
cmd <- paste0('mv ', topDir, '*.maf ', mutdir)
system(cmd)



