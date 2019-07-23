# generate patient report
rmarkdown::render(input = 'OMPARE.Rmd', 
                  params = list(topDir = 'data/Test/',
                                fusion_method = 'arriba',
                                set_title = 'Test Report'))

