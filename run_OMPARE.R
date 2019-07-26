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
