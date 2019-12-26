# Step2: generate patient report

# parameter info:
# fusion_method can be either arriba, star, both or not specified
# fusion_method and set_title is completely optional 
# topDir is mandatory

# some examples on how to run OMPARE
# 1. fusion_method is arriba
rmarkdown::render(input = 'OMPARE.Rmd', 
                  params = list(topDir = 'data/PNOC008/',
                                fusion_method = 'arriba',
                                set_title = 'Patient Report'))


# 2. fusion_method is star
rmarkdown::render(input = 'OMPARE.Rmd', 
                  params = list(topDir = 'data/PNOC008/',
                                fusion_method = 'star',
                                set_title = 'Patient Report'))


# 3. fusion_method is both
rmarkdown::render(input = 'OMPARE.Rmd', 
                  params = list(topDir = 'data/PNOC008/',
                                fusion_method = 'both',
                                set_title = 'Patient Report'))


# 4. fusion_method not specified (uses both)
rmarkdown::render(input = 'OMPARE.Rmd', 
                  params = list(topDir = 'data/PNOC008/',
                                set_title = 'Patient Report'))

# PNOC008-04
# only consensus calls
rmarkdown::render(input = 'OMPARE.Rmd', 
                  params = list(topDir = 'data/PNOC008-04/',
                                fusion_method = 'arriba',
                                set_title = 'PNOC008-04 Patient Report',
                                snv_consensus = TRUE),
                  output_file = 'PNOC008_04_consensus.html')

# all four callers
rmarkdown::render(input = 'OMPARE.Rmd', 
                  params = list(topDir = 'data/PNOC008-04/',
                                fusion_method = 'arriba',
                                set_title = 'PNOC008-04 Patient Report',
                                snv_consensus = FALSE),
                  output_file = 'PNOC008_04.html')
